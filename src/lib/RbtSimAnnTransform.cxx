/***********************************************************************
 * The rDock program was developed from 1998 - 2006 by the software team
 * at RiboTargets (subsequently Vernalis (R&D) Ltd).
 * In 2006, the software was licensed to the University of York for
 * maintenance and distribution.
 * In 2012, Vernalis and the University of York agreed to release the
 * program as Open Source software.
 * This version is licensed under GNU-LGPL version 3.0 with support from
 * the University of Barcelona.
 * http://rdock.sourceforge.net/
 ***********************************************************************/

// Simple, non-adaptive simulated annealing protocol
#include <iomanip>
using std::setw;

#include "RbtBaseSF.h"
#include "RbtChrom.h"
#include "RbtSFRequest.h"
#include "RbtSimAnnTransform.h"
#include "RbtWorkSpace.h"

// Simple class to keep track of Monte Carlo sampling statistics
RbtMCStats::RbtMCStats() { Init(0.0); }

void RbtMCStats::Init(RbtDouble score) {
    _min = _max = _initial = _final = score;
    InitBlock(score);
}

void RbtMCStats::InitBlock(RbtDouble score) {
    _blockMin = _blockMax = _blockInitial = _blockFinal = score;
    _total = 0.0;
    _total2 = 0.0;
    _steps = 0;
    _accepted = 0;
}

void RbtMCStats::Accumulate(RbtDouble score, RbtBool bAccepted) {
    _steps++;
    if (bAccepted) _accepted++;
    _total += score;
    _total2 += score * score;
    _blockMin = std::min(_blockMin, score);
    _blockMax = std::max(_blockMax, score);
    _blockFinal = _final = score;
    _min = std::min(_min, score);
    _max = std::max(_max, score);
}

RbtDouble RbtMCStats::Mean() const { return _total / _steps; }
RbtDouble RbtMCStats::Variance() const { return _total2 / _steps - pow(Mean(), 2); }
RbtDouble RbtMCStats::AccRate() const { return float(_accepted) / float(_steps); }

// Static data member for class type
RbtString RbtSimAnnTransform::_CT("RbtSimAnnTransform");
// Parameter names
RbtString RbtSimAnnTransform::_START_T("START_T");
RbtString RbtSimAnnTransform::_FINAL_T("FINAL_T");
RbtString RbtSimAnnTransform::_BLOCK_LENGTH("BLOCK_LENGTH");
RbtString RbtSimAnnTransform::_SCALE_CHROM_LENGTH("SCALE_CHROM_LENGTH");
RbtString RbtSimAnnTransform::_NUM_BLOCKS("NUM_BLOCKS");
RbtString RbtSimAnnTransform::_STEP_SIZE("STEP_SIZE");
RbtString RbtSimAnnTransform::_MIN_ACC_RATE("MIN_ACC_RATE");
RbtString RbtSimAnnTransform::_PARTITION_DIST("PARTITION_DIST");
RbtString RbtSimAnnTransform::_PARTITION_FREQ("PARTITION_FREQ");
RbtString RbtSimAnnTransform::_HISTORY_FREQ("HISTORY_FREQ");

const RbtSimAnnTransform::Config
    RbtSimAnnTransform::DEFAULT_CONFIG{};  // Empty initializer to fall back to default values

RbtSimAnnTransform::RbtSimAnnTransform(const RbtString& strName, const Config& config):
    RbtBaseBiMolTransform(_CT, strName),
    m_rand(Rbt::GetRbtRand()),
    config{config}  // For now simply copy the config structure.
{
    m_spStats = RbtMCStatsPtr(new RbtMCStats());
    DEBUG_ERR(_CT << " parameterised constructor" << endl);
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSimAnnTransform::~RbtSimAnnTransform() {
    DEBUG_ERR(_CT << " destructor" << endl);
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

////////////////////////////////////////
// Protected methods
///////////////////
void RbtSimAnnTransform::SetupReceptor() {}

void RbtSimAnnTransform::SetupLigand() {}

void RbtSimAnnTransform::SetupTransform() {
    // Construct the overall chromosome for the system
    m_chrom.SetNull();
    m_lastGoodVector.clear();
    m_minVector.clear();
    RbtWorkSpace* pWorkSpace = GetWorkSpace();
    if (pWorkSpace) {
        m_chrom = new RbtChrom(pWorkSpace->GetModels());
        RbtInt chromLength = m_chrom->GetLength();
        m_lastGoodVector.reserve(chromLength);
        m_minVector.reserve(chromLength);
    }
}

////////////////////////////////////////
// Private methods
///////////////////
// Pure virtual in RbtBaseTransform
// Actually apply the transform
void RbtSimAnnTransform::Execute() {
    // Get the current scoring function from the workspace
    RbtWorkSpace* pWorkSpace = GetWorkSpace();
    if (pWorkSpace == NULL)  // Return if this transform is not registered
        return;
    RbtBaseSF* pSF = pWorkSpace->GetSF();
    if (pSF == NULL)  // Return if workspace does not have a scoring function
        return;
    RbtInt iTrace = GetTrace();

    pWorkSpace->ClearPopulation();

    RbtInt block_length = config.block_length;
    if (config.scale_chromosome_length) {
        RbtInt chromLength = m_chrom->GetLength();
        block_length *= chromLength;
    }
    if (iTrace > 0) cout << _CT << ": Block length = " << block_length << endl;

    // Send the partitioning request separately based on the current partition distance
    // If partDist is zero, the partitioning automatically gets removed
    m_spPartReq = RbtRequestPtr(new RbtSFPartitionRequest(config.partition_distance));
    pSF->HandleRequest(m_spPartReq);

    // Update the chromosome to match the current model coords
    // and keep track of the vector corresponding to minimum score
    m_chrom->SyncFromModel();
    m_minVector.clear();
    m_chrom->GetVector(m_minVector);
    RbtDouble score = pSF->Score();
    m_spStats->Init(score);

    cout.precision(3);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout.setf(ios_base::right, ios_base::adjustfield);

    if (iTrace > 0) {
        cout << _CT << ": Initial score = " << score << endl;
        cout << endl
             << endl
             << setw(5) << "BLOCK" << setw(10) << "TEMP" << setw(10) << "ACC.RATE" << setw(10) << "STEP" << setw(10)
             << "INITIAL" << setw(10) << "FINAL" << setw(10) << "MEAN" << setw(10) << "S.DEV" << setw(10) << "MIN"
             << setw(10) << "MAX" << endl;
    }

    RbtDouble initial_temp = config.initial_temp;
    RbtDouble step_size = config.step_size;

    // DM 15 Feb 1999 - don't initialise the Monte Carlo stats each block
    // if we are doing a constant temperature run
    RbtBool bInitBlock = (config.initial_temp != config.final_temp);
    RbtDouble tFac =
        (config.num_blocks > 1) ? pow(config.final_temp / config.initial_temp, 1.0 / (config.num_blocks - 1)) : 1.0;

    for (RbtInt iBlock = 1; iBlock <= config.num_blocks; iBlock++, initial_temp *= tFac) {
        if (bInitBlock) {
            m_spStats->InitBlock(pSF->Score());
        }
        MC(initial_temp, block_length, step_size);
        if (iTrace > 0) {
            cout << setw(5) << iBlock << setw(10) << initial_temp << setw(10) << m_spStats->AccRate() << setw(10)
                 << step_size << setw(10) << m_spStats->_blockInitial << setw(10) << m_spStats->_blockFinal << setw(10)
                 << m_spStats->Mean() << setw(10) << sqrt(m_spStats->Variance()) << setw(10) << m_spStats->_blockMin
                 << setw(10) << m_spStats->_blockMax << endl;
        }

        // Halve the maximum step sizes for all enabled modes
        // if the acceptance rate is less than the threshold
        if (m_spStats->AccRate() < config.min_accuracy_rate) {
            step_size *= 0.5;
            // Reinitialise the stats (only need to do it here if bInitBlock is false
            // otherwise it will be done at the beginning of the next block)
            if (!bInitBlock) {
                m_spStats->InitBlock(pSF->Score());
            }
        }
    }
    // Update the model coords with the minimum score chromosome
    m_chrom->SetVector(m_minVector);
    m_chrom->SyncToModel();
    RbtRequestPtr spClearPartReq(new RbtSFPartitionRequest(0.0));
    pSF->HandleRequest(spClearPartReq);  // Clear any partitioning
    if (iTrace > 0) {
        cout << endl << _CT << ": Final score = " << pSF->Score() << endl;
    }
}

void RbtSimAnnTransform::MC(RbtDouble t, RbtInt blockLen, RbtDouble stepSize) {
    RbtBaseSF* pSF = GetWorkSpace()->GetSF();
    RbtInt iTrace = GetTrace();
    RbtDouble score = pSF->Score();

    // Keep a record of the last good chromosome vector, for fast revert following a failed Metropolic test
    m_lastGoodVector.clear();
    m_chrom->GetVector(m_lastGoodVector);
    // Main loop over number of MC steps
    for (RbtInt iStep = 1; iStep <= blockLen; iStep++) {
        m_chrom->Mutate(stepSize);
        m_chrom->SyncToModel();
        RbtDouble newScore = pSF->Score();
        RbtDouble delta = newScore - score;
        RbtBool bMetrop = ((delta < 0.0) || (exp(-1000.0 * delta / (8.314 * t)) > m_rand.GetRandom01()));
        // PASSED
        if (bMetrop) {
            score = newScore;
            // Update the last good vector
            m_lastGoodVector.clear();
            m_chrom->GetVector(m_lastGoodVector);
            // Update the minimum score vector
            if (score < m_spStats->_min) {
                m_minVector.clear();
                m_chrom->GetVector(m_minVector);
            }
        }
        // FAILED
        else {
            // revert to old chromosome
            // No need to SyncToModel as this will be done after the next mutation
            m_chrom->SetVector(m_lastGoodVector);
        }
        // Gather the statistics
        m_spStats->Accumulate(score, bMetrop);
        // Render to the history file if appropriate (true = with component scores)
        if ((config.history_frequency > 0) && (iStep % config.history_frequency) == 0) {
            GetWorkSpace()->SaveHistory(true);
        }

        // Update the interaction lists if appropriate
        // We update the lists every nth accepted trial (as rejected trials don't change the coords)
        if ((config.partition_frequency > 0) && ((m_spStats->_accepted % config.partition_frequency) == 0)) {
            pSF->HandleRequest(m_spPartReq);
            RbtDouble oldScore = score;
            score = pSF->Score();
            if (fabs(score - oldScore) > 0.001 && (iTrace > 1)) {
                cout << "** WARNING - Interaction lists updated, change in score = " << score - oldScore << endl;
            }
        }
    }
    // Ensure the model is synchronised with the chromosome on exit
    // This would not be the case if the last Metropolis test were failed
    m_chrom->SyncToModel();
}
