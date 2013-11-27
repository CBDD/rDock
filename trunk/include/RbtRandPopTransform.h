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

//Generates a randomised population for the GA
#ifndef _RBTRANDPOPTRANSFORM_H_
#define _RBTRANDPOPTRANSFORM_H_

#include "RbtBaseBiMolTransform.h"
#include "RbtChromElement.h"

class RbtRandPopTransform : public RbtBaseBiMolTransform
{
 public:
  static RbtString _CT;
  static RbtString _POP_SIZE;
  static RbtString _SCALE_CHROM_LENGTH;
    
  ////////////////////////////////////////
  //Constructors/destructors
  RbtRandPopTransform(const RbtString& strName = "RANDPOP");
  virtual ~RbtRandPopTransform();

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
  RbtRandPopTransform(const RbtRandPopTransform&);//Copy constructor disabled by default
  RbtRandPopTransform& operator=(const RbtRandPopTransform&);//Copy assignment disabled by default

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
typedef SmartPtr<RbtRandPopTransform> RbtRandPopTransformPtr;//Smart pointer

#endif //_RBTRANDPOPTRANSFORM_H_
