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

#include "RbtSFAgg.h"
#include "RbtWorkSpace.h"
#include "RbtModel.h"

//Static data member for class type
RbtString RbtSFAgg::_CT("RbtSFAgg");

////////////////////////////////////////
//Constructors/destructors
RbtSFAgg::RbtSFAgg(const RbtString& strName) : RbtBaseSF(_CT,strName),m_nNonHLigandAtoms(0) {
#ifdef _DEBUG
	cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSFAgg::~RbtSFAgg() {
	//Delete all our children
#ifdef _DEBUG
	cout << _CT << "::~" << _CT << "(): Deleting child scoring functions of " << GetName() << endl;
#endif //_DEBUG
	//We need to iterate using a while loop because each deletion will reduce the size of m_sf,
	//hence conventional iterators would become invalid
	while (m_sf.size() > 0) {
		RbtBaseSF* pSF = m_sf.back();
		//Assertion: parent of child is this object
		Assert<RbtAssert>(!SFAGG_CHECK||pSF->m_parent==this);
#ifdef _DEBUG
		cout << "Deleting " << pSF->GetName() << " from " << GetName() << endl;
#endif //_DEBUG
		delete pSF;
	}
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}


////////////////////////////////////////
//Public methods
////////////////

//Returns all child component scores as a string-variant map
//Key = fully qualified component name, value = weighted score
//(for saving in a Model's data fields)
void RbtSFAgg::ScoreMap(RbtStringVariantMap& scoreMap) const {
  if (isEnabled()) {
    //First populate all the child entries in the scoring function map
    for (RbtBaseSFListConstIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
        (*iter)->ScoreMap(scoreMap);
    }
    //Now we can pick up the raw score for this aggregate from the map entry
    //without calculating it directly,
    //as the child terms will have incremented the total as part of their ScoreMap
    //implementation
    RbtString name = GetFullName();
    RbtDouble rs = scoreMap[name];
    
    //Cascade to the parent of this aggregate
    AddToParentMapEntry(scoreMap, rs);

    //7 Feb 2005 (DM, Enspiral Discovery) - normalise the raw aggregate score by
    //the number of heavy (non-H) ligand atoms
    if (m_nNonHLigandAtoms > 0) {
      RbtDouble norm_s = rs / m_nNonHLigandAtoms;
      scoreMap[name+".norm"] = norm_s;
      //Only record the number of heavy atoms for the root aggregate (SCORE)
      //root = has no parent
      if (!GetParentSF()) {
        scoreMap[name+".heavy"] = m_nNonHLigandAtoms;
      }
    }
  }
}

//Aggregate handling methods
void RbtSFAgg::Add(RbtBaseSF* pSF) throw (RbtError) {
	//By first orphaning the scoring function to be added,
	//we handle attempts to readd existing children automatically,
	pSF->Orphan();
	pSF->m_parent = this;
#ifdef _DEBUG
		cout << _CT << "::Add(): Adding " << pSF->GetName() << " to " << GetName() << endl;
#endif //_DEBUG
	m_sf.push_back(pSF);
}

void RbtSFAgg::Remove(RbtBaseSF* pSF) throw (RbtError) {
	RbtBaseSFListIter iter = std::find(m_sf.begin(),m_sf.end(),pSF);
	if (iter == m_sf.end()) {
		throw RbtBadArgument(_WHERE_,"Remove(): pSF is not a member of this aggregate");
	}
	else {
		//Assertion: parent of child is this object
		Assert<RbtAssert>(!SFAGG_CHECK||pSF->m_parent==this);
#ifdef _DEBUG
		cout << _CT << "::Remove(): Removing " << pSF->GetName() << " from " << GetName() << endl;
#endif //_DEBUG
		m_sf.erase(iter);
		pSF->m_parent = NULL;//Nullify the parent pointer of the child that has been removed 
	}
}

RbtBool RbtSFAgg::isAgg() const {return true;}
RbtUInt RbtSFAgg::GetNumSF() const {return m_sf.size();}
RbtBaseSF* RbtSFAgg::GetSF(RbtUInt iSF) const throw (RbtError)
{
	if (iSF >= m_sf.size()) {
		throw RbtBadArgument(_WHERE_,"GetSF(): iSF out of range");
	}
	else {
		return m_sf[iSF];
	}
}

//WorkSpace handling methods
//Register scoring function with a workspace
//Aggregate version registers all children, AND itself (new behaviour, 7 Feb 2005 (DM))
void RbtSFAgg::Register(RbtWorkSpace* pWorkSpace) {
	//7 Feb 2005 (DM, Enspiral Discovery) - register the aggregate as well
	//as we need the number of heavy atoms in the ligand for calculating
	//the normalised scores
	RbtBaseObject::Register(pWorkSpace);
#ifdef _DEBUG
		cout << _CT << "::Register(): Registering child scoring functions of " << GetName() << endl;
#endif //_DEBUG
	for (RbtBaseSFListIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
		(*iter)->Register(pWorkSpace);
	}
}

//Unregister with a workspace
//Aggregate version unregisters all children, AND itself (new behaviour, 7 Feb 2005 (DM))
void RbtSFAgg::Unregister() {
#ifdef _DEBUG
		cout << _CT << "::Unregister(): Unregistering child scoring functions of " << GetName() << endl;
#endif //_DEBUG
	for (RbtBaseSFListIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
		(*iter)->Unregister();
	}
	RbtBaseObject::Unregister();//7 Feb 2005 (DM) - unregister the agg, for the same reasons as above
}

//Override RbtObserver pure virtual
//Notify observer that subject has changed
//7 Feb 2005 (DM, Enspiral Discovery) - keep track of number of non-H ligand atoms
//so can report aggregate scores normalised by #heavy atoms
void RbtSFAgg::Update(RbtSubject* theChangedSubject) {
  RbtWorkSpace* pWorkSpace = GetWorkSpace();
  if (theChangedSubject == pWorkSpace) {
    //Check if ligand has been updated (model #1)
    if (pWorkSpace->GetNumModels() >= 2) {
      RbtModelPtr spLigand = pWorkSpace->GetModel(1);
      if (spLigand.Ptr()) {
        m_nNonHLigandAtoms = Rbt::GetNumAtoms(spLigand->GetAtomList(),std::not1(Rbt::isAtomicNo_eq(1)));
      }
      else {
        m_nNonHLigandAtoms = 0;
      }
    }
  }
}

//Request Handling method
//Aggregates handle the request themselves first, then cascade to all children
void RbtSFAgg::HandleRequest(RbtRequestPtr spRequest) {
	RbtBaseObject::HandleRequest(spRequest);
	for (RbtBaseSFListIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
		(*iter)->HandleRequest(spRequest);
	}
}

//Virtual function for dumping scoring function details to an output stream
//Called by operator <<
void RbtSFAgg::Print(ostream& s) const {
	//First print the parameters for this aggregate
	RbtBaseObject::Print(s);
	//Now call Print() for each of the children
	for (RbtBaseSFListConstIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
		(*iter)->Print(s);
	}
}

////////////////////////////////////////
//Private methods
/////////////////
//Raw score for an aggregate is the sum of the weighted scores of its children
RbtDouble RbtSFAgg::RawScore() const {
	RbtDouble score(0.0);
	for (RbtBaseSFListConstIter iter = m_sf.begin(); iter != m_sf.end(); iter++) {
		score += (*iter)->Score();
	}
	return score;
}

