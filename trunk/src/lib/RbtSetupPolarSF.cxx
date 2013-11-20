/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtSetupPolarSF.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include <iomanip>
using std::setw;
#include "RbtSetupPolarSF.h"

//Static data members
RbtString RbtSetupPolarSF::_CT("RbtSetupPolarSF");
RbtString RbtSetupPolarSF::_RADIUS("RADIUS");
RbtString RbtSetupPolarSF::_NORM("NORM");
RbtString RbtSetupPolarSF::_POWER("POWER");
RbtString RbtSetupPolarSF::_CHGFACTOR("CHGFACTOR");
RbtString RbtSetupPolarSF::_GUANFACTOR("GUANFACTOR");
	
//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtSetupPolarSF::RbtSetupPolarSF(const RbtString& strName)
  : RbtBaseSF(_CT,strName)
{
  //Add parameters
  AddParameter(_RADIUS,5.0);
  AddParameter(_NORM,25);
  AddParameter(_POWER,0.5);
  AddParameter(_CHGFACTOR,0.5);
  AddParameter(_GUANFACTOR,1.0);
  Disable();
#ifdef _DEBUG
  cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtSetupPolarSF::~RbtSetupPolarSF() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtSetupPolarSF::SetupReceptor() {
  if (GetReceptor().Null())
    return;
  RbtAtomList atomList = GetReceptor()->GetAtomList();
  RbtAtomList heavyAtomList = Rbt::GetAtomList(atomList,std::not1(Rbt::isAtomicNo_eq(1)));
  RbtInt traceTriggerLevel = 1;
  SetupAtomList(atomList,heavyAtomList,traceTriggerLevel);
}

void RbtSetupPolarSF::SetupLigand() {
  if (GetLigand().Null())
    return;
  RbtAtomList atomList = GetLigand()->GetAtomList();
  RbtAtomList emptyList;
  RbtInt traceTriggerLevel = 2;
  SetupAtomList(atomList, emptyList, traceTriggerLevel);
}

void RbtSetupPolarSF::SetupSolvent() {
  RbtModelList solventList = GetSolvent();
  if (solventList.empty())
    return;
  //Concatenate all solvent atoms into a single list
  RbtAtomList solventAtomList;
  for (RbtModelListConstIter iter = solventList.begin(); iter != solventList.end(); ++iter) {
  	RbtAtomList atomList = (*iter)->GetAtomList();
   	std::copy(atomList.begin(),atomList.end(),std::back_inserter(solventAtomList));
  }
  RbtAtomList emptyList;
  RbtInt traceTriggerLevel = 1;
  SetupAtomList(solventAtomList, emptyList, traceTriggerLevel);
}

//Reusable method for setting up an arbitrary atom list
//If neigbourList is empty, then do not include the neighbour density factor
void RbtSetupPolarSF::SetupAtomList(RbtAtomList& atomList,
									const RbtAtomList& neighbourList,
									RbtInt traceTriggerLevel) {
  if (atomList.empty())
    return;
  
  RbtInt iTrace = GetTrace();
  RbtBool bCalcNeighbourDensity = !neighbourList.empty();
  RbtDouble radius = GetParameter(_RADIUS);
  RbtDouble norm = GetParameter(_NORM);
  RbtDouble power = GetParameter(_POWER);
  RbtDouble chgFactor = GetParameter(_CHGFACTOR);
  RbtDouble guanFactor = GetParameter(_GUANFACTOR);

  Rbt::isAtomHBondAcceptor bIsHBA;
  Rbt::isAtomHBondDonor bIsHBD;
  Rbt::isAtomGuanidiniumCarbon bIsGuan;
  Rbt::isAtomLipophilic bIsLipo;

  std::ios_base::fmtflags oldflags = cout.flags();
  cout.precision(3);
  cout.setf(ios_base::fixed,ios_base::floatfield);
  cout.setf(ios_base::right,ios_base::adjustfield);
  if (iTrace > traceTriggerLevel) {
    cout << endl << _CT << ": Atom\t\t, f(neighb), f(charge), isLipo,  isHBA,  isHBD" << endl;
  }
  
  for (RbtAtomListIter iter = atomList.begin(); iter != atomList.end(); iter++) {
  	RbtDouble fNeighb = 1.0;
  	if (bCalcNeighbourDensity) {
	    //Neighbour density function
	    RbtInt nNeighb = Rbt::GetNumAtoms(neighbourList,Rbt::isAtomInsideSphere((*iter)->GetCoords(),radius));
	    nNeighb--;//Don't count the atom itself!
	    fNeighb = pow(nNeighb/norm,power);
  	}
    //Charge factor (always >= 1, sign reflects whether charge is positive or negative)
    RbtDouble charge = (*iter)->GetGroupCharge();
    //DM 18 Jan 2001 - need to check for neutral HBond acceptors and set sign to -1
    RbtDouble sign = ( (charge < 0.0) || (bIsHBA(*iter)) ) ? -1.0 : 1.0;
    charge = sign*(1.0+(fabs(charge)*chgFactor));
    if (bIsGuan(*iter)) {
      charge *= guanFactor;//Adjustable weight for "IONIC" interactions with guanidinium carbons
    }
// XB apply reweighting factor for receptor's polar atoms:
    if (traceTriggerLevel == 1) {
      RbtDouble wxb = (*iter)->GetReweight();
      (*iter)->SetUser1Value(fNeighb*charge*wxb);
     // cout << "Atom: " << (*iter)->GetAtomName() << " Polar_User1Value: " << (*iter)->GetUser1Value() << " weigth: " << wxb << endl;
    }else{
      (*iter)->SetUser1Value(fNeighb*charge);
    }
// XB END MODIFICATIONS
    (*iter)->SetUser1Flag(bIsLipo(*iter));
    if (iTrace > traceTriggerLevel) {
   //   cout << _CT << ": " << (*iter)->GetFullAtomName() << "\t," << setw(10) << fNeighb << "," << setw(10) << charge
//	   << "," << setw(7) << (*iter)->GetUser1Flag() << "," << setw(7) << bIsHBA(*iter) << "," << setw(7) << bIsHBD(*iter) << endl;
    }
  }
  cout.flags(oldflags);
}

void RbtSetupPolarSF::SetupScore() {
  //No further setup required
}

RbtDouble RbtSetupPolarSF::RawScore() const {
  return 0.0;
}

