/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPopulation.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Manages a population of genomes.
//Main public method (GAstep) performs one iteration of a GA.
#ifndef _RBTPOPULATION_H_
#define _RBTPOPULATION_H_

#include "RbtGenome.h"
#include "RbtError.h"

class RbtBaseSF; //forward definition

class RbtPopulation {
public:
  static RbtString _CT;
  //Sole constructor to create a randomised genome population of a fixed size.
  //pChr is the seed chromosome to clone to create each genome.
  //size is the population size to create.
  //pSF is the scoring function used to rank the genomes.
  //Notes:
  //1) Each genome is randomised by calling the Randomise() method
  //of the underlying chromosome.
  //2) Scoring function scores are calculated.
  //3) Genomes are sorted by score.
  //4) Roulette wheel fitness values are calculated.
  //5) Model coords are updated to match the fittest chromosome
  //An RbtBadArgument error is thrown if size is <=0, or if pChr or pSF is null.
  RbtPopulation(RbtChromElement* pChr, RbtInt size, RbtBaseSF* pSF) throw (RbtError);
  virtual ~RbtPopulation();
  
  //Gets the maximum size of the population as defined in the constructor.
  RbtInt GetMaxSize() const {return m_size;}
  //Gets the actual size of the population (may be < GetMaxSize())
  RbtInt GetActualSize() const {return m_pop.size();}
  //Gets the best genome in the sorted population (element zero).
  RbtGenomePtr Best() const;
  //Gets the average raw score across entire population
  RbtDouble GetScoreMean() const {return m_scoreMean;}
  //Gets the raw score variance across entire population
  RbtDouble GetScoreVariance() const {return m_scoreVariance;}
  //Gets the vector of genomes in the population.
  const RbtGenomeList& GetGenomeList() const {return m_pop;}
  
  //Gets the scoring function used for ranking genomes
  RbtBaseSF* GetSF() const {return m_pSF;}
  //Sets the scoring function used for ranking genomes
  //Forces a recalculation of the genome scores, and the roulette wheel fitness values
  //SetSF should be called whenever the scoring function parameters have changed
  //e.g. in between GA stages.
  //An RbtBadArgument error is thrown if pSF is null.
  //Model coords are updated to match the fittest chromosome
  void SetSF(RbtBaseSF* pSF) throw (RbtError);
  
  //Main method for performing a GA iteration
  void GAstep(RbtInt nReplicates,//Number of new genomes to create in the iteration
                RbtDouble relStepSize,//Relative step size for chromosome mutations
                RbtDouble equalityThreshold,//Equality threshold for genomes 
                RbtDouble pcross,//Probability of crossover
                RbtBool xovermut,//if true, perform cauchy mutation following crossover
                RbtBool cmutate //true=cauchy mutations, false=regular mutations
                ) throw (RbtError);
  RbtGenomePtr RouletteWheelSelect() const;
  
  void Print(ostream&) const;
  friend ostream& operator<<(ostream& , const RbtPopulation &);
  
private:
  //Merges the new individuals created into the main population
  //Duplicate genomes are removed (based on equality of chromosome elements, not scores)
  void MergeNewPop(RbtGenomeList& newPop, RbtDouble equalityThreshold);
  void EvaluateRWFitness();
  RbtPopulation(const RbtPopulation&);//Disable
  RbtPopulation& operator=(const RbtPopulation &);//Disable
  
  RbtGenomeList m_pop;//The population of genomes
  RbtInt m_size;//The maximum size of the population
  RbtDouble m_c;//Sigma Truncation Multiplier
  RbtBaseSF* m_pSF;//The scoring function
  RbtRand& m_rand;//reference to the singleton random number generator
  RbtDouble m_scoreMean;//the average raw score across all genomes
  RbtDouble m_scoreVariance;//the variance of raw scores across all genomes
};

typedef SmartPtr<RbtPopulation> RbtPopulationPtr;

#endif  //_RBTPOPULATION_H_
 
