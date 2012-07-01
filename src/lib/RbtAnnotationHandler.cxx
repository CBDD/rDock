/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

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
