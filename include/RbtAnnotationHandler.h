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

//Base implementation class for managing an annotation list. Scoring function
//classes wishing to store annotations for later retrieval should derive from
//this class.

#ifndef _RBTANNOTATIONHANDLER_H_
#define _RBTANNOTATIONHANDLER_H_

#include "RbtAnnotation.h"

class RbtAnnotationHandler
{
 public:
  static RbtString _ANNOTATION_FIELD;

  ////////////////////////////////////////
  //Constructors/destructors
  virtual ~RbtAnnotationHandler(); //Default destructor  
  
  ////////////////////////////////////////
  //Public methods
  ////////////////
  RbtBool isAnnotationEnabled() const;
  //Get a const ref to the annotation list (read-only)
  const RbtAnnotationList& GetAnnotationList() const;
  RbtInt GetNumAnnotations() const;    
  void RenderAnnotationList(const RbtString& strName, RbtStringList& retVal) const;
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  RbtAnnotationHandler();
  void AddAnnotation(RbtAnnotationPtr spAnnotation) const;
  void ClearAnnotationList() const;
  void EnableAnnotations(RbtBool bEnabled) const;
  
 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtAnnotationHandler(const RbtAnnotationHandler&);//Copy constructor disabled by default  
  RbtAnnotationHandler& operator=(const RbtAnnotationHandler&);//Copy assignment disabled by default
  
  
 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////
  
  
 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  mutable RbtBool m_bEnabled;
  mutable RbtAnnotationList m_annotationList;
};

#endif //_RBTANNOTATIONHANDLER_H_
