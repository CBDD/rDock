/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtBaseTransform.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Abstract base class for all classes which manipulate molecular coords

#ifndef _RBTBASETRANSFORM_H_
#define _RBTBASETRANSFORM_H_

#include "RbtConfig.h"
#include "RbtBaseObject.h"

class RbtTransformAgg;//forward declaration

class RbtBaseTransform : public RbtBaseObject
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names

  ////////////////////////////////////////
  //Constructors/destructors
  virtual ~RbtBaseTransform();
  
  //Give aggregates access to RbtBaseSF private methods/data
  friend class RbtTransformAgg;

  ////////////////////////////////////////
  //Public methods
  ////////////////
  //Fully qualified name, prefixed by all ancestors
  RbtString GetFullName() const;

  //Main public method - actually apply the transform
  //Not virtual. Base class method checks if transform is enabled,
  //sends SFRequests to reconfigure the scoring function
  //then calls private virtual method Execute() to apply the transform
  void Go();
  
  //Aggregate handling methods
  virtual void Add(RbtBaseTransform*) throw (RbtError);
  virtual void Remove(RbtBaseTransform*) throw (RbtError);
  virtual RbtBool isAgg() const;
  virtual RbtUInt GetNumTransforms() const;
  virtual RbtBaseTransform* GetTransform(RbtUInt) const throw (RbtError);
  void Orphan();//Force removal from the parent aggregate
  RbtBaseTransform* GetParentTransform() const;
	
  //Scoring function request handling
  //Transforms can store up a list of requests to send to the workspace
  //scoring function each time Go() is executed
  void AddSFRequest(RbtRequestPtr);
  void ClearSFRequests();
  void SendSFRequests();
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  RbtBaseTransform(const RbtString& strClass, const RbtString& strName);

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtBaseTransform();
  RbtBaseTransform(const RbtBaseTransform&);//Copy constructor disabled by default
  RbtBaseTransform& operator=(const RbtBaseTransform&);//Copy assignment disabled by default
  
  //PURE VIRTUAL - subclasses should override. Applies the transform
  virtual void Execute() = 0;
  
 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtBaseTransform* m_parent;
  RbtRequestList m_SFRequests;
};

//Useful typedefs
typedef vector<RbtBaseTransform*> RbtBaseTransformList;//Vector of smart pointers
typedef RbtBaseTransformList::iterator RbtBaseTransformListIter;
typedef RbtBaseTransformList::const_iterator RbtBaseTransformListConstIter;

#endif //_RBTBASETRANSFORM_H_
