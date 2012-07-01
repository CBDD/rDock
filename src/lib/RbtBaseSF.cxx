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

#include "RbtBaseSF.h"
#include "RbtSFRequest.h"

//Static data members
RbtString RbtBaseSF::_CT("RbtBaseSF");
RbtString RbtBaseSF::_WEIGHT("WEIGHT");
RbtString RbtBaseSF::_RANGE("RANGE");
RbtString RbtBaseSF::_SYSTEM_SF("SCORE.SYSTEM");
RbtString RbtBaseSF::_INTRA_SF("SCORE.INTRA");

////////////////////////////////////////
//Constructors/destructors
RbtBaseSF::RbtBaseSF(const RbtString& strClass, const RbtString& strName)
  : RbtBaseObject(strClass,strName), m_parent(NULL), m_weight(1.0), m_range(10.0) {
#ifdef _DEBUG
  cout << _CT << " parameterised constructor for " << strClass << endl;
#endif //_DEBUG
  //Add parameters
  AddParameter(_WEIGHT,m_weight);
  AddParameter(_RANGE,m_range);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

//Dummy default constructor for virtual base subclasses
//Should never get called
RbtBaseSF::RbtBaseSF() {
#ifdef _DEBUG
  cout << "WARNING: " << _CT << " default constructor" << endl;
#endif //_DEBUG
  //_RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtBaseSF::~RbtBaseSF() {
  Orphan();//Remove object from parent aggregate
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}


////////////////////////////////////////
//Public methods
////////////////

//Fully qualified name, prefixed by all ancestors (e.g. SCORE.INTER.HBOND)
RbtString RbtBaseSF::GetFullName() const {
  //If we have a parent, prefix our short name with our parents full name,
  //else full name is just our short name
  return (m_parent) ? m_parent->GetFullName() + "." + GetName() : GetName();
}

RbtDouble RbtBaseSF::GetWeight() const {return m_weight;}
void RbtBaseSF::SetWeight(RbtDouble weight) {SetParameter(_WEIGHT,weight);}

RbtDouble RbtBaseSF::GetRange() const {return m_range;}
void RbtBaseSF::SetRange(RbtDouble range) {SetParameter(_RANGE,range);}

//Returns weighted score if scoring function is enabled, else returns zero
RbtDouble RbtBaseSF::Score() const {
  return isEnabled() ? GetWeight()*RawScore() : 0.0;
}

//Returns all child component scores as a string-variant map
//Key = fully qualified component name, value = weighted score
//(for saving in a Model's data fields)
void RbtBaseSF::ScoreMap(RbtStringVariantMap& scoreMap) const {
  if (isEnabled()) {
    //DM 17 Jan 2006.
    //New approach:
    //1) We record the raw, unweighted score for this term
    //in the map
    RbtDouble rs = RawScore();
    RbtString name = GetFullName();
    scoreMap[name] = rs;
    //2) We add the weighted score to the parent aggregate
    //entry. This gives us the opportunity to override
    //ScoreMap in order to divert scores away from their
    //natural parent entry. e.g. SCORE.INTER.VDW could
    //record intra-receptor and intra-solvent contributions
    //under SCORE.SYSTEM.VDW
    AddToParentMapEntry(scoreMap, rs);
  }
}

//Helper method for ScoreMap
void RbtBaseSF::AddToParentMapEntry(RbtStringVariantMap& scoreMap, RbtDouble rs) const {
    if (m_parent) {
        RbtDouble w = GetWeight();
        RbtDouble s = w * rs;
        RbtString parentName = m_parent->GetFullName();
        RbtDouble parentScore = scoreMap[parentName];
        parentScore += s;
        scoreMap[parentName] = parentScore;
    }
}
		
//Aggregate handling (virtual) methods
//Base class throws an InvalidRequest error

void RbtBaseSF::Add(RbtBaseSF*) throw (RbtError) {
  throw RbtInvalidRequest(_WHERE_,"Add() invalid for non-aggregate scoring functions");
}
void RbtBaseSF::Remove(RbtBaseSF*) throw (RbtError) {
  throw RbtInvalidRequest(_WHERE_,"Remove() invalid for non-aggregate scoring functions");
}

RbtBaseSF* RbtBaseSF::GetSF(RbtUInt iSF) const throw (RbtError) {
  throw RbtInvalidRequest(_WHERE_,"GetSF() invalid for non-aggregate scoring functions");
}
RbtBool RbtBaseSF::isAgg() const {return false;}
RbtUInt RbtBaseSF::GetNumSF() const {return 0;}

//Aggregate handling (concrete) methods
RbtBaseSF* RbtBaseSF::GetParentSF() const {return m_parent;}
//Force removal from the parent aggregate
void RbtBaseSF::Orphan() {
  if (m_parent) {
#ifdef _DEBUG
    cout << _CT << "::Orphan(): Removing " << GetName() << " from " << m_parent->GetName() << endl;
#endif //_DEBUG
    m_parent->Remove(this);
  }
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtBaseSF::ParameterUpdated(const RbtString& strName) {
  if (strName == _WEIGHT) {
    m_weight = GetParameter(_WEIGHT);
  }
  else if (strName == _RANGE) {
    m_range = GetParameter(_RANGE);
  }
  else {
    RbtBaseObject::ParameterUpdated(strName);
  }
}

