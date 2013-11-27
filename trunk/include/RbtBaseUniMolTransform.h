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
