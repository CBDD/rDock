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

#include "RbtAnnotationHandler.h"

//Static data members
RbtString RbtAnnotationHandler::_ANNOTATION_FIELD("ANNOTATIONS");

////////////////////////////////////////
//Constructors/destructors
RbtAnnotationHandler::RbtAnnotationHandler() : m_bEnabled(false) {
  _RBTOBJECTCOUNTER_CONSTR_("RbtAnnotationHandler");
}

RbtAnnotationHandler::~RbtAnnotationHandler() {
  _RBTOBJECTCOUNTER_DESTR_("RbtAnnotationHandler");
}

////////////////////////////////////////
//Public methods
////////////////
RbtBool RbtAnnotationHandler::isAnnotationEnabled() const {
  return m_bEnabled;
}

const RbtAnnotationList& RbtAnnotationHandler::GetAnnotationList() const {
  return m_annotationList;
}

RbtInt RbtAnnotationHandler::GetNumAnnotations() const {
  return m_annotationList.size();
}

void RbtAnnotationHandler::RenderAnnotationList(const RbtString& strName, RbtStringList& retVal) const {
  //It is callers responsibility to clear retVal before calling RenderAnnotationList
  //Allows for appending to existing string list
  for (RbtAnnotationListConstIter iter = m_annotationList.begin(); iter != m_annotationList.end(); iter++) {
    retVal.push_back(strName + "," + (*iter)->Render());
  }
}

////////////////////////////////////////
//Protected methods
///////////////////
void RbtAnnotationHandler::AddAnnotation(RbtAnnotationPtr spAnnotation) const {
  if (m_bEnabled) {
    m_annotationList.push_back(spAnnotation);
  }
}

void RbtAnnotationHandler::ClearAnnotationList() const {
  m_annotationList.clear();
}

void RbtAnnotationHandler::EnableAnnotations(RbtBool bEnabled) const {
  m_bEnabled = bEnabled;
}
