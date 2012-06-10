/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtAnnotationHandler.h#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
