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

#include "RbtGATransform.h"

#include <iomanip>

#include "RbtPopulation.h"
#include "RbtSFRequest.h"
#include "RbtWorkSpace.h"
using std::setw;

RbtString RbtGATransform::_CT("RbtGATransform");
RbtString RbtGATransform::_NEW_FRACTION("NEW_FRACTION");
RbtString RbtGATransform::_PCROSSOVER("PCROSSOVER");
RbtString RbtGATransform::_XOVERMUT("XOVERMUT");
RbtString RbtGATransform::_CMUTATE("CMUTATE");
RbtString RbtGATransform::_STEP_SIZE("STEP_SIZE");
RbtString RbtGATransform::_EQUALITY_THRESHOLD("EQUALITY_THRESHOLD");
RbtString RbtGATransform::_NCYCLES("NCYCLES");
RbtString RbtGATransform::_NCONVERGENCE("NCONVERGENCE");
RbtString RbtGATransform::_HISTORY_FREQ("HISTORY_FREQ");

const RbtGATransform::Config RbtGATransform::DEFAULT_CONFIG{};  // Empty initializer to fall back to default values

RbtGATransform::RbtGATransform(const RbtString& strName, const Config& config):
    RbtBaseBiMolTransform(_CT, strName),
    m_rand(Rbt::GetRbtRand()),
    config{config} {
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtGATransform::~RbtGATransform() { _RBTOBJECTCOUNTER_DESTR_(_CT); }

void RbtGATransform::SetupReceptor() {}

void RbtGATransform::SetupLigand() {}

void RbtGATransform::SetupTransform() {}

void RbtGATransform::Execute() {
    RbtWorkSpace* pWorkSpace = GetWorkSpace();
    if (pWorkSpace == NULL) {
        return;
    }
    RbtBaseSF* pSF = pWorkSpace->GetSF();
    if (pSF == NULL) {
        return;
    }
    RbtPopulationPtr pop = pWorkSpace->GetPopulation();
    if (pop.Null() || (pop->GetMaxSize() < 1)) {
        return;
    }
    // Remove any partitioning from the scoring function
    // Not appropriate for a GA
    pSF->HandleRequest(new RbtSFPartitionRequest(0.0));
    // This forces the population to rescore all the individuals in case
    // the scoring function has changed
    pop->SetSF(pSF);

    RbtInt popsize = pop->GetMaxSize();
    RbtInt nrepl = config.population_size_fraction_as_new_individuals_per_cycle * popsize;
    RbtBool bHistory = config.history_frequency > 0;
    RbtInt iTrace = GetTrace();

    RbtDouble bestScore = pop->Best()->GetScore();
    // Number of consecutive cycles with no improvement in best score
    RbtInt iConvergence = 0;

    if (iTrace > 0) {
        cout.precision(3);
        cout.setf(ios_base::fixed, ios_base::floatfield);
        cout.setf(ios_base::right, ios_base::adjustfield);
        cout << endl
             << setw(5) << "CYCLE" << setw(5) << "CONV" << setw(10) << "BEST" << setw(10) << "MEAN" << setw(10)
             << "VAR" << endl;

        cout << endl
             << setw(5) << "Init" << setw(5) << "-" << setw(10) << bestScore << setw(10) << pop->GetScoreMean()
             << setw(10) << pop->GetScoreVariance() << endl;
    }

    for (RbtInt iCycle = 0; (iCycle < config.max_cycles) && (iConvergence < config.num_convergence_cycles); ++iCycle) {
        if (bHistory && ((iCycle % config.history_frequency) == 0)) {
            pop->Best()->GetChrom()->SyncToModel();
            pWorkSpace->SaveHistory(true);
        }
        pop->GAstep(
            nrepl,
            config.relative_step_size,
            config.equality_threshold,
            config.crossover_probability,
            config.cauchy_mutation_after_crossover,
            config.use_cauchy_distribution_for_mutations
        );
        RbtDouble score = pop->Best()->GetScore();
        if (score > bestScore) {
            bestScore = score;
            iConvergence = 0;
        } else {
            iConvergence++;
        }
        if (iTrace > 0) {
            cout << setw(5) << iCycle << setw(5) << iConvergence << setw(10) << score << setw(10)
                 << pop->GetScoreMean() << setw(10) << pop->GetScoreVariance() << endl;
        }
    }
    pop->Best()->GetChrom()->SyncToModel();
    RbtInt ri = GetReceptor()->GetCurrentCoords();
    GetLigand()->SetDataValue("RI", ri);
}
