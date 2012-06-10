/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtBaseIntraSF.h#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Virtual base class for intramolecular scoring functions. Overrides Update()

#ifndef _RBTBASEINTRASF_H_
#define _RBTBASEINTRASF_H_

#include "RbtBaseSF.h"
#include "RbtModel.h"

class RbtBaseIntraSF : public virtual RbtBaseSF
{
 public:
  //Class type string
  static RbtString _CT;

  ////////////////////////////////////////
  //Constructors/destructors
  virtual ~RbtBaseIntraSF();
  
  ////////////////////////////////////////
  //Public methods
  ////////////////
  
  RbtModelPtr GetLigand() const;
  
  //Override RbtObserver pure virtual
  //Notify observer that subject has changed
  virtual void Update(RbtSubject* theChangedSubject);
  
  //Override RbtBaseSF::ScoreMap to provide additional raw descriptors
  virtual void ScoreMap(RbtStringVariantMap& scoreMap) const;
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  RbtBaseIntraSF();
  
  //PURE VIRTUAL - Derived classes must override
  virtual void SetupScore() = 0;//Called by Update when model has changed
  
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
  RbtModelPtr m_spLigand;
  //26 Mar 2003 (DM) Remember the raw score for the initial ligand conformation
  //This becomes the zero point for all subsequent score reporting
  //i.e. all intramolecular scores are reported relative to the initial score
  RbtDouble m_zero;
};

#endif //_RBTBASEINTRASF_H_
