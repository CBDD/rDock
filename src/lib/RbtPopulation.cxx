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

#include "RbtPopulation.h"
#include "RbtDebug.h"
#include "RbtDockingError.h"
#include <algorithm>

RbtString RbtPopulation::_CT("RbtPopulation"); 

RbtPopulation::RbtPopulation(RbtChromElement* pChr, RbtInt size, RbtBaseSF* pSF)
    throw (RbtError)
        : m_size(size), m_c(2.0), m_pSF(pSF), m_rand(Rbt::GetRbtRand()),
        m_scoreMean(0.0), m_scoreVariance(0.0)
{
  if (pChr == NULL) {
    throw RbtBadArgument(_WHERE_, "Null chromosome element passed to RbtPopulation constructor");
  }
  else if (size <= 0) {
    throw RbtBadArgument(_WHERE_, "Population size must be positive (non-zero)");
  }
  //Create a random population
  m_pop.reserve(m_size);
  for (RbtInt i = 0; i < m_size; ++i) {
    //The RbtGenome constructor clones the chromosome to create an independent copy
    RbtGenomePtr genome = new RbtGenome(pChr);
    genome->GetChrom()->Randomise();
    m_pop.push_back(genome);
  }
  //Calculate the scores and evaluate roulette wheel fitness
  SetSF(m_pSF);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtPopulation::~RbtPopulation() {
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

//Sets the scoring function used for ranking genomes
void RbtPopulation::SetSF(RbtBaseSF* pSF) throw (RbtError) {
  if (pSF == NULL) {
    throw RbtBadArgument(_WHERE_, "Null scoring function passed to SetSF");
  }
  m_pSF = pSF;
  for (RbtGenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    (*iter)->SetScore(m_pSF);
  }
  std::stable_sort(m_pop.begin(), m_pop.end(), Rbt::GenomeCmp_Score());
  EvaluateRWFitness();
}

void RbtPopulation::GAstep(RbtInt nReplicates, RbtDouble relStepSize,
                            RbtDouble equalityThreshold, RbtDouble pcross,
                            RbtBool xovermut, RbtBool cmutate) throw (RbtError)
{
  if (nReplicates <= 0) {
    throw RbtBadArgument(_WHERE_, "nReplicates must be positive (non-zero)");
  }
  RbtGenomeList newPop;
  newPop.reserve(nReplicates);
  for (RbtInt i = 0 ; i < nReplicates / 2 ; i++) {
    RbtGenomePtr mother = RouletteWheelSelect();
    RbtGenomePtr father = RouletteWheelSelect();
    //Check that mother and father are not the same genome
    //The check is on the pointers, not that the chromosomes are near-equal
    //If we repeatedly get the same genomes selected, this must mean
    //the population lacks diversity
    RbtInt j = 0;
    while (father == mother) {
      father = RouletteWheelSelect();
      if (j > 100)  
        throw RbtDockingError(_WHERE_, 
				            "Population failure - not enough diversity");
      j++;
    }
    RbtGenomePtr child1 = new RbtGenome(*mother);
    RbtGenomePtr child2 = new RbtGenome(*father);
    //Crossover
    if (m_rand.GetRandom01() < pcross) {
      Rbt::Crossover(father->GetChrom(), mother->GetChrom(), 
                     child1->GetChrom(), child2->GetChrom());
      //Cauchy mutation following crossover
      if (xovermut) {
        child1->GetChrom()->CauchyMutate(0.0, relStepSize);
        child2->GetChrom()->CauchyMutate(0.0, relStepSize);
      }
    }
    //Mutation
    else {
      //Cauchy mutation
      if (cmutate) {
        child1->GetChrom()->CauchyMutate(0.0, relStepSize);
        child2->GetChrom()->CauchyMutate(0.0, relStepSize);
      }
      //Regular mutation
      else {
        child1->GetChrom()->Mutate(relStepSize);
        child2->GetChrom()->Mutate(relStepSize);
      }
    }
    newPop.push_back(child1);
    newPop.push_back(child2);
  }
  //check if one more is needed (odd nReplicates).
  if (nReplicates % 2) {
    RbtGenomePtr mother = RouletteWheelSelect();
    RbtGenomePtr child = new RbtGenome(*mother);
    child->GetChrom()->CauchyMutate(0.0, relStepSize);
    newPop.push_back(child);
  }
  MergeNewPop(newPop, equalityThreshold);
  EvaluateRWFitness();
}

RbtGenomePtr RbtPopulation::Best() const {
  return m_pop.empty() ? RbtGenomePtr() : m_pop.front();
}

void RbtPopulation::Print(ostream& s) const {
  s << m_pop.size() << endl;
  for (RbtGenomeListConstIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    (*iter)->Print(s);
  }
}

ostream& operator<<(ostream& s, const RbtPopulation& p) {
  p.Print(s);
  return s;
}


RbtGenomePtr RbtPopulation::RouletteWheelSelect() const {
  RbtDouble cutoff = m_rand.GetRandom01();
  RbtInt size = m_pop.size();
  RbtInt lower = 0;
  RbtInt upper = size - 1;
  while (upper >= lower) {
    RbtInt i = lower + (upper - lower) / 2;
    if (m_pop[i]->GetRWFitness() > cutoff)
      upper = i - 1;
    else
      lower = i + 1;
  } 
  //make sure lower is a number between 0 and size - 1
  lower = std::min(size-1, lower);
  lower = std::max(0, lower);
  return (m_pop[lower]);
}

void RbtPopulation::MergeNewPop(RbtGenomeList& newPop, RbtDouble equalityThreshold) {
  //Assume newPop needs scoring and sorting
  for (RbtGenomeListIter iter = newPop.begin(); iter != newPop.end(); ++iter) {
    (*iter)->SetScore(m_pSF);
  }
  std::stable_sort(newPop.begin(), newPop.end(), Rbt::GenomeCmp_Score());
     
  RbtGenomeList mergedPop;
  mergedPop.reserve(m_pop.size() + newPop.size());
  //Merge pops by score
  std::merge(m_pop.begin(), m_pop.end(),
                newPop.begin(), newPop.end(), 
                std::back_inserter(mergedPop), Rbt::GenomeCmp_Score());
  //Remove neighbouring duplicates by equality of chromosome element values
  RbtGenomeListIter end = std::unique(mergedPop.begin(), mergedPop.end(),
                                        Rbt::isGenome_eq(equalityThreshold));
  mergedPop.erase(end, mergedPop.end());
  m_pop.clear();
  end = (mergedPop.size() > m_size) ? (mergedPop.begin() + m_size) : mergedPop.end();
  std::copy(mergedPop.begin(), end, back_inserter(m_pop));
}

void RbtPopulation::EvaluateRWFitness() {
  //Determine mean and variance of true scores
  RbtDouble sum(0.0);
  RbtDouble sumSq(0.0);
  for (RbtGenomeListConstIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    RbtDouble score = (*iter)->GetScore();
    sum += score;
    sumSq += score*score;
  }
  RbtDouble popSize = m_pop.size();
  m_scoreMean = sum / popSize;
  m_scoreVariance = (sumSq / popSize) - (m_scoreMean * m_scoreMean);
  RbtDouble sigma = sqrt(m_scoreVariance);
  // calculate scaled fitness values using sigma truncation
  // Goldberg page 124
  RbtDouble offset = m_scoreMean - m_c * sigma;
  RbtDouble partialSum = 0.0;
  //First set the unnormalised fitness values
  for (RbtGenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    partialSum = (*iter)->SetRWFitness(offset, partialSum);
  }
  //Now normalise so that fitness values run from 0 to 1, ready for
  //roulette wheel selection
  for (RbtGenomeListIter iter = m_pop.begin(); iter != m_pop.end(); ++iter) {
    (*iter)->NormaliseRWFitness(partialSum);
  }
}

