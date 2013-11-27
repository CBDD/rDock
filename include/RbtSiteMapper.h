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

//Base class for docking site mapper function classes
//Designed in a similar fashion to scoring functions, in that mapper objects
//are attached to a workspace and can respond to changes in the defined
//receptor and ligand models
//Sub classes should implement operator() which returns a list of cavity smart
//pointers

#ifndef _RBTSITEMAPPER_H_
#define _RBTSITEMAPPER_H_

#include "RbtConfig.h"
#include "RbtBaseObject.h"
#include "RbtCavity.h"
#include "RbtModel.h"

class RbtSiteMapper : public RbtBaseObject
{
  public:
  //Class type string
  static RbtString _CT;

  ////////////////////////////////////////
  //Constructors/destructors
  virtual ~RbtSiteMapper();
  ////////////////////////////////////////
  //Public methods
  ////////////////
  RbtModelPtr GetReceptor() const {return m_spReceptor;}

  //PURE VIRTUAL - subclasses must override
  //NB - subclasses should also override RbtObserver::Update pure virtual
  virtual RbtCavityList operator() () = 0;

  //Override RbtObserver pure virtual
  //Notify observer that subject has changed
  virtual void Update(RbtSubject* theChangedSubject);

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  RbtSiteMapper(const RbtString& strClass, const RbtString& strName);

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtSiteMapper();//default constructor disabled
  RbtSiteMapper(const RbtSiteMapper&);//Copy constructor disabled by default      
  RbtSiteMapper& operator=(const RbtSiteMapper&);//Copy assignment disabled by default

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////
      
 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtModelPtr m_spReceptor;
};
  
//Useful typedefs
typedef SmartPtr<RbtSiteMapper> RbtSiteMapperPtr;//Smart pointer

#endif //_RBTSITEMAPPER_H_
