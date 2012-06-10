/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtBaseUniMolTransform.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Base class for unimolecular transforms. Overrides Update()

#ifndef _RBTBASEUNIMOLTRANSFORM_H_
#define _RBTBASEUNIMOLTRANSFORM_H_

#include "RbtBaseTransform.h"
#include "RbtModel.h"

class RbtBaseUniMolTransform : public RbtBaseTransform
{
 public:
  //Class type string
  static RbtString _CT;
  ////////////////////////////////////////
  //Constructors/destructors
  virtual ~RbtBaseUniMolTransform();
  
  ////////////////////////////////////////
  //Public methods
  ////////////////
  
  RbtModelPtr GetLigand() const;
  
  //Override RbtObserver pure virtual
  //Notify observer that subject has changed
  virtual void Update(RbtSubject* theChangedSubject);
  
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  RbtBaseUniMolTransform(const RbtString& strClass, const RbtString& strName);
  
  //PURE VIRTUAL - Derived classes must override
  virtual void SetupTransform() = 0;//Called by Update when model has changed
  
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
};

#endif //_RBTBASEUNIMOLTRANSFORM_H_
