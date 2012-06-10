/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtAnnotationHandler.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
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
