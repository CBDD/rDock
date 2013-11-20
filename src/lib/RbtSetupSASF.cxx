/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtSetupSASF.cxx#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
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

