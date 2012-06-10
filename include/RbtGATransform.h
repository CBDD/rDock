/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtGATransform.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
