/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSimplexTransform.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Simplex Minimiser
#ifndef _RBTSIMPLEXTRANSFORM_H_
#define _RBTSIMPLEXTRANSFORM_H_

#include "RbtBaseBiMolTransform.h"
#include "RbtChromElement.h"

class RbtSimplexTransform : public RbtBaseBiMolTransform
{
 public:
  //Static data member for class type
  static RbtString _CT;
    // Parameter names
  static RbtString _MAX_CALLS;
  static RbtString _NCYCLES;
  static RbtString _STOPPING_STEP_LENGTH;
  static RbtString _PARTITION_DIST;
  static RbtString _STEP_SIZE;
  //Stop once score improves by less than convergence value
  //between cycles
  static RbtString _CONVERGENCE;
   
  ////////////////////////////////////////
  //Constructors/destructors
  RbtSimplexTransform(const RbtString& strName = "SIMPLEX");
  virtual ~RbtSimplexTransform();

  ////////////////////////////////////////
  //Public methods
  ////////////////	
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  virtual void SetupTransform();//Called by Update when either model has changed
  virtual void SetupReceptor();  // Called by Update when receptor is changed
  virtual void SetupLigand();    // Called by Update when ligand is changed
  virtual void SetupSolvent();    // Called by Update when solvent is changed
  virtual void Execute();

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtSimplexTransform(const RbtSimplexTransform&);//Copy constructor disabled by default
  RbtSimplexTransform& operator=(const RbtSimplexTransform&);//Copy assignment disabled by default
 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtChromElementPtr m_chrom;
};

//Useful typedefs
typedef SmartPtr<RbtSimplexTransform> RbtSimplexTransformPtr;//Smart pointer
	
#endif //_RBTSIMPLEXTRANSFORM_H_
