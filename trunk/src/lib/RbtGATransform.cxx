/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

#include "RbtGATransform.h"
#include "RbtWorkSpace.h"
#include "RbtPopulation.h"
#include "RbtSFRequest.h"
#include <iomanip>
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

RbtGATransform::RbtGATransform(const RbtString& strName) :
                               RbtBaseBiMolTransform(_CT,strName),
                               m_rand(Rbt::GetRbtRand())
{
  AddParameter(_NEW_FRACTION, 0.5);
  AddParameter(_PCROSSOVER, 0.4);
  AddParameter(_XOVERMUT, true);
  AddParameter(_CMUTATE, false);
  AddParameter(_STEP_SIZE, 1.0);
  AddParameter(_EQUALITY_THRESHOLD, 0.1);
  AddParameter(_NCYCLES, 100);
  AddParameter(_NCONVERGENCE, 6);
  AddParameter(_HISTORY_FREQ, 0);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtGATransform::~RbtGATransform() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtGATransform::SetupReceptor(){}

void RbtGATransform::SetupLigand(){}

void RbtGATransform::SetupTransform(){}

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
  //Remove any partitioning from the scoring function
  //Not appropriate for a GA
  pSF->HandleRequest(new RbtSFPartitionRequest(0.0));
  //This forces the population to rescore all the individuals in case
  //the scoring function has changed
  pop->SetSF(pSF);
  
  RbtDouble newFraction = GetParameter(_NEW_FRACTION);
  RbtDouble pcross = GetParameter(_PCROSSOVER);
  RbtBool xovermut = GetParameter(_XOVERMUT);
  RbtBool cmutate = GetParameter(_CMUTATE);
  RbtDouble relStepSize = GetParameter(_STEP_SIZE);
  RbtDouble equalityThreshold = GetParameter(_EQUALITY_THRESHOLD);
  RbtInt nCycles = GetParameter(_NCYCLES);
  RbtInt nConvergence = GetParameter(_NCONVERGENCE);
  RbtInt nHisFreq = GetParameter(_HISTORY_FREQ);
  
  RbtInt popsize = pop->GetMaxSize();
  RbtInt nrepl = newFraction * popsize;
  RbtBool bHistory = nHisFreq > 0;
  RbtInt iTrace = GetTrace();

  RbtDouble bestScore = pop->Best()->GetScore();
   //Number of consecutive cycles with no improvement in best score
  RbtInt iConvergence = 0;
 
  if (iTrace > 0) {
    cout.precision(3);
    cout.setf(ios_base::fixed,ios_base::floatfield);
    cout.setf(ios_base::right,ios_base::adjustfield);
    cout << endl
         << setw(5) << "CYCLE"
         << setw(5) << "CONV"
         << setw(10) << "BEST"
         << setw(10) << "MEAN"
         << setw(10) << "VAR"
         << endl;
         
    cout << endl
         << setw(5)  << "Init"
         << setw(5)  << "-"
         << setw(10) << bestScore
         << setw(10) << pop->GetScoreMean()
         << setw(10) << pop->GetScoreVariance()
         << endl;
  }

  for (RbtInt iCycle = 0;
                (iCycle < nCycles) && (iConvergence < nConvergence);
                ++iCycle) {
    if (bHistory && ( (iCycle % nHisFreq) == 0)){
        pop->Best()->GetChrom()->SyncToModel();
        pWorkSpace->SaveHistory(true); 
    }
    pop->GAstep(nrepl, relStepSize, equalityThreshold, pcross, xovermut, cmutate);
    RbtDouble score = pop->Best()->GetScore();
    if (score > bestScore) {
        bestScore = score;
        iConvergence = 0;
    }
    else {
        iConvergence++;
    }
    if (iTrace > 0) {
        cout << setw(5)  << iCycle
             << setw(5)  << iConvergence
             << setw(10) << score
             << setw(10) << pop->GetScoreMean()
             << setw(10) << pop->GetScoreVariance()
             << endl;
    }
  }
  pop->Best()->GetChrom()->SyncToModel();
  RbtInt ri = GetReceptor()->GetCurrentCoords();
  GetLigand()->SetDataValue("RI",ri);
}
