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

//Evolves an existing population using a GA
#ifndef _RBTGATRANSFORM_H_
#define _RBTGATRANSFORM_H_

#include "RbtBaseBiMolTransform.h"
#include "RbtRand.h"
#include "RbtPopulation.h"

class RbtGATransform : public RbtBaseBiMolTransform {
 public:
  static RbtString _CT;
  //New individuals to create each cycle, as fraction of population size
  static RbtString _NEW_FRACTION;
  //Probability of crossover (1-probability of mutation)
  static RbtString _PCROSSOVER;
  //If true, perform Cauchy mutation after each crossover
  static RbtString _XOVERMUT;
  //If true, mutations are from Cauchy distribution; if false, from rect. distribution
  static RbtString _CMUTATE;
  //Relative step size for mutations (relative to absolute step sizes defined
  //for each chromosome element)
  static RbtString _STEP_SIZE;
  //Two genomes are considered equal if the maximum relative difference
  //between chromosome elements is less than _EQUALITY_THRESHOLD
  static RbtString _EQUALITY_THRESHOLD;
  //Maximum number of cycles
  static RbtString _NCYCLES;
  //Terminate if the best score does not improve over _NCONVERGENCE
  //consecutive cycles
  static RbtString _NCONVERGENCE;
  //Output the best pose every _HISTORY_FREQ cycles.
  static RbtString _HISTORY_FREQ;


    ////////////////////////////////////////
    //Constructors/destructors
    ////////////////////////////////////////
  RbtGATransform(const RbtString& strName = "GAGENRW");
  virtual ~RbtGATransform();

   
 protected:
    ////////////////////////////////////////
    //Protected methods
    ///////////////////
  virtual void SetupReceptor();  //Called by Update when receptor is changed
  virtual void SetupLigand();   //Called by Update when ligand is changed
  virtual void SetupTransform();//Called by Update when either model has changed
  virtual void Execute();

 private:
    ////////////////////////////////////////
    //Private methods
    /////////////////
  RbtGATransform(const RbtGATransform&);//Copy constructor disabled by default
  RbtGATransform& operator=(const RbtGATransform&);//Copy assignment disabled by default

 private:
  RbtRand& m_rand;
}; 

#endif //_RBTGATRANSFORM_H_
