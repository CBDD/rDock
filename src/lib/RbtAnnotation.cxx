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

#include "RbtAnnotation.h"
#include "RbtPseudoAtom.h"

////////////////////////////////////////
//Constructors/destructors
RbtAnnotation::RbtAnnotation(const RbtAtom* pAtom1, const RbtAtom* pAtom2, RbtDouble dist, RbtDouble score) :
  m_pAtom1(pAtom1), m_pAtom2(pAtom2), m_dist(dist), m_score(score) {
}

RbtAnnotation::~RbtAnnotation() {
}
  
////////////////////////////////////////
//Public methods
////////////////
//Get
const RbtAtom* RbtAnnotation::GetAtom1Ptr() const {return m_pAtom1;}
const RbtAtom* RbtAnnotation::GetAtom2Ptr() const {return m_pAtom2;}
RbtDouble RbtAnnotation::GetDistance() const {return m_dist;}
RbtDouble RbtAnnotation::GetScore() const {return m_score;}

//Get the fully qualified (FQ) residue name for atom 2 (target atom)
RbtString RbtAnnotation::GetFQResName() const {
  return m_pAtom2->GetSegmentName() + ":" + m_pAtom2->GetSubunitName() + "_" + m_pAtom2->GetSubunitId() + ":";
}

//Set
void RbtAnnotation::SetAtom1Ptr(const RbtAtom* pAt1) {m_pAtom1 = pAt1;}
void RbtAnnotation::SetAtom2Ptr(const RbtAtom* pAt2) {m_pAtom2 = pAt2;}
void RbtAnnotation::SetDistance(RbtDouble d) {m_dist = d;}
void RbtAnnotation::SetScore(RbtDouble s) {m_score = s;}

//Render annotation into string in rDock Viewer format
//Note: string does not contain the first field (annotation name)
//Full string required for SD file is:
//strName + "," + Render();
RbtString RbtAnnotation::Render() const {
  ostrstream ostr;
  //Check if either atom is a pseudoatom
  //If so, arbitrarily write the first non-bridgehead atom ID
  const RbtPseudoAtom* pseudo1 = dynamic_cast<const RbtPseudoAtom*>(m_pAtom1);
  const RbtPseudoAtom* pseudo2 = dynamic_cast<const RbtPseudoAtom*>(m_pAtom2);
  if (pseudo1) {
    RbtAtomList tmpList = pseudo1->GetAtomList();
    RbtAtomList al1 = Rbt::GetAtomList(tmpList,std::not1(Rbt::isAtomBridgehead()));
    //Check for rare case of multifused rings, where all constituent ring atoms are bridgeheads
    //In this case, just use the first atom in the ring list
    RbtAtomPtr spAt1 = (al1.empty()) ? tmpList.front() : al1.front();
    ostr << spAt1->GetAtomId() << ",";
  }
  else {
    ostr << m_pAtom1->GetAtomId() << ",";
  }
  if (pseudo2) {
    RbtAtomList tmpList = pseudo2->GetAtomList();
    RbtAtomList al2 = Rbt::GetAtomList(tmpList,std::not1(Rbt::isAtomBridgehead()));
    //Check for rare case of multifused rings, where all constituent ring atoms are bridgeheads
    //In this case, just use the first atom in the ring list
    RbtAtomPtr spAt2 = (al2.empty()) ? tmpList.front() : al2.front();
    ostr << spAt2->GetAtomId() << ",";
  }
  else {
    ostr << m_pAtom2->GetAtomId() << ",";
  }
  ostr.setf(ios_base::fixed,ios_base::floatfield);
  ostr.precision(2);
  ostr << m_dist << "," << m_score << ends;
  RbtString retVal(ostr.str());
  delete ostr.str();
  return retVal;
}

//Operators
//Special meaning of operator+ for accumulating annotations by residue
void RbtAnnotation::operator+=(const RbtAnnotation& ann) {
  //Check FQ residue names match, if not then don't accumulate
  if (GetFQResName() != ann.GetFQResName()) {
    return;
  }
  //atom 1 and distance relate to the closest contact between the ligand and the
  //target residue, so update them here if necessary
  if (ann.GetDistance() < m_dist) {
    m_pAtom1 = ann.GetAtom1Ptr();
    m_dist = ann.GetDistance();
  }
  //Score represents the total score between the ligand and the target residue
  m_score += ann.GetScore();
  return;
}
