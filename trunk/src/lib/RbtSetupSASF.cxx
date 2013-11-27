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

#include "RbtSetupSASF.h"

RbtString RbtSetupSASF::_CT("RbtSetupSASF");

RbtSetupSASF::RbtSetupSASF(const RbtString& strName) : RbtBaseSF(_CT, strName)
{
#ifdef _DEBUG
	cout << _CT << "parameterized constructor      <--------" << endl;
#endif	//_DEBUG
	Disable();
	_RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSetupSASF::~RbtSetupSASF()
{
#ifdef _DEBUG
	cout << _CT << " destructor                    <--------" << endl;
#endif	//_DEBUG
	_RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtSetupSASF::SetupReceptor()
{
#ifdef _DEBUG
	cout << _CT << " SetupReceptor                 <--------" << endl;
#endif	//_DEBUG
	// get receptor atoms 
	//theReceptorList = Rbt::GetAtomList(GetReceptor()->GetAtomList(),std::not1(Rbt::isAtomicNo_eq(1)));
	theReceptorList = GetReceptor()->GetAtomList();
#ifdef _DEBUG
	cout << _CT << "::SetupReceptor(): #ATOMS = " << theReceptorList.size() << endl;
#endif	//_DEBUG
}

void RbtSetupSASF::SetupScore()
{
#ifdef _DEBUG
	cout << _CT << " SetupScore                    <--------" << endl;
#endif	//_DEBUG
}

void RbtSetupSASF::SetupLigand()
{
#ifdef _DEBUG
	cout << _CT << " SetupLigand                   <--------" << endl;
#endif	//_DEBUG
	theLigandList = GetLigand()->GetAtomList();
#ifdef _DEBUG
	cout << _CT << "::SetupLigand(): #ATOMS = " << theLigandList.size() << endl;
#endif	//_DEBUG
}

RbtDouble RbtSetupSASF::RawScore() const { return 0.0;}

