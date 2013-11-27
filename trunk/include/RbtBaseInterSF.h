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

//Virtual base class for intermolecular scoring functions. Overrides Update()

#ifndef _RBTBASEINTERSF_H_
#define _RBTBASEINTERSF_H_

#include "RbtBaseSF.h"
#include "RbtModel.h"

class RbtBaseInterSF : public virtual RbtBaseSF
{
 public:
  //Class type string
  static RbtString _CT;

  ////////////////////////////////////////
  //Constructors/destructors
  virtual ~RbtBaseInterSF();
  
  ////////////////////////////////////////
  //Public methods
  ////////////////
  
  RbtModelPtr GetReceptor() const;
  RbtModelPtr GetLigand() const;
  RbtModelList GetSolvent() const;
  
  //Override RbtObserver pure virtual
  //Notify observer that subject has changed
  virtual void Update(RbtSubject* theChangedSubject);
  
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  RbtBaseInterSF();
  
  //PURE VIRTUAL - Derived classes must override
  virtual void SetupReceptor() = 0;//Called by Update when receptor is changed
  virtual void SetupLigand() = 0;//Called by Update when ligand is changed
  virtual void SetupSolvent() {};//Called by Update when solvent is changed
  virtual void SetupScore() = 0;//Called by Update when either model has changed
  
 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  
 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////
  
  
 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtModelPtr m_spReceptor;
  RbtModelPtr m_spLigand;
  RbtModelList m_solventList;
};

#endif //_RBTBASEINTERSF_H_
