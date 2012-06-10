/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSiteMapper.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
