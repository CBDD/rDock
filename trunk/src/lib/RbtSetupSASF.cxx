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

