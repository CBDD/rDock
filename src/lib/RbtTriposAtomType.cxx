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

#include "RbtTriposAtomType.h"
#include "RbtAtom.h"
#include "RbtBond.h"

////////////////////////////////////////
//Constructors/destructors

RbtTriposAtomType::RbtTriposAtomType() {
  SetupTypeInfo();
}

RbtTriposAtomType::~RbtTriposAtomType() {

}

////////////////////////////////////////
//Public methods
////////////////

//Returns Tripos atom type for given atom
//If useExtendedTypes is true, pseudo types for extended carbons and polar Hs are used
RbtTriposAtomType::eType RbtTriposAtomType::operator() (RbtAtom* pAtom,RbtBool useExtendedTypes) const {
  RbtInt nAtomicNo = pAtom->GetAtomicNo();
  RbtAtom::eHybridState eState = pAtom->GetHybridState();
  RbtInt nImplH = pAtom->GetNumImplicitHydrogens();
  RbtBondMap bondMap;
  Rbt::isBondAmide bIsAmide;
  Rbt::isAtomHBondDonor bIsDonor;
  Rbt::isAtomGuanidiniumCarbon bIsGuan;

  switch (nAtomicNo) {

  case 1://H
    if (useExtendedTypes && bIsDonor(pAtom))
      return H_P;
    else
      return H;

  case 3://Li
    return Li;

  case 6://C
    //DM 19 Jun 2003 - C.cat is the Tripos type for guanidinium carbon
    if (bIsGuan(pAtom))
      return C_cat;

    switch (eState) {

    case RbtAtom::SP:
      if (useExtendedTypes)
	switch (nImplH) {
	case 1: return C_1_H1;
	default:return C_1;
	}
      else
	return C_1;

    case RbtAtom::SP2:
      if (useExtendedTypes)
	switch (nImplH) {
	case 1: return C_2_H1;
	case 2: return C_2_H2;
	default:return C_2;
	}
      else
	return C_2;

    case RbtAtom::SP3:
      if (useExtendedTypes)
	switch (nImplH) {
	case 1: return C_3_H1;
	case 2: return C_3_H2;
	case 3: return C_3_H3;
	default:return C_3;
	}
      else
	return C_3;

    case RbtAtom::AROM:
      if (useExtendedTypes)
	switch (nImplH) {
	case 1: return C_ar_H1;
	default:return C_ar;
	}
      else
	return C_ar;
      
    default:
      return C_3;
    }

  case 7://N
    switch (eState) {
    case RbtAtom::SP:
      return N_1;
    case RbtAtom::SP2:
      return N_2;
    case RbtAtom::SP3:
      if (pAtom->GetNumBonds() == 4)
	return N_4;
      else
	return N_3;
    case RbtAtom::AROM:
      return N_ar;
    case RbtAtom::TRI:
      bondMap = pAtom->GetBondMap();
      for (RbtBondMapConstIter iter = bondMap.begin(); iter != bondMap.end(); iter++) {
	if (bIsAmide((*iter).first)) { // passing the bond
	  return N_am;
	}
      }
      return N_pl3;
    default:
      return N_3;
    }

  case 8://O
    //DM 19 Jun 2003 - O.co2 is the Tripos type for neg charged oxygen (e.g. in acid groups)
    if (pAtom->GetGroupCharge() < 0.0)
      return O_co2;
    
    switch (eState) {
    case RbtAtom::SP2:
      return O_2;
    case RbtAtom::SP3:
    case RbtAtom::TRI:
      //DM 7 Nov 2003 - trap nitro oxygens here (one is O_TRI, the other is O_SP2,
      //both should be typed as O_2 for vdW purposes)
      return (pAtom->GetNumBonds() == 1) ? O_2 : O_3;
    default:
      return O_3;
    }
    
  case 9://F
    return F;

  case 11://Na
    return Na;

  case 13://Al
    return Al;

  case 14://Si
    return Si;

  case 15://P
    return P_3;

  case 16://S
    switch (eState) {
    case RbtAtom::SP2:
      switch (pAtom->GetCoordinationNumber("O_SP2")) {
      case 2:
	return S_o2;
      case 1:
	return S_o;
      default:
	return S_2;
      }
    case RbtAtom::SP3:
      return S_3;
    default:
      return S_3;
    }

  case 17://Cl
    return Cl;
    
  case 19://K
    return K;

  case 20://Ca
    return Ca;

  case 35://Br
    return Br;

  case 53://I
    return I;

  default:
    return  UNDEFINED;
  }
}

//Converts Tripos type to string
RbtString RbtTriposAtomType::Type2Str(RbtTriposAtomType::eType aType) const {
  return m_typeInfo[aType].name;
}
// Get hybridisation from Tripos type
int RbtTriposAtomType::Type2Hybrid(RbtTriposAtomType::eType aType) const {
  return m_typeInfo[aType].hybrid;
}
// Get atomic number from Tripos type
RbtInt RbtTriposAtomType::Type2AtomicNo(RbtTriposAtomType::eType aType) const {
  return m_typeInfo[aType].atomicNo;
}

//Converts string to Tripos type
RbtTriposAtomType::eType RbtTriposAtomType::Str2Type(const RbtString& strType) const {
  for (RbtTriposInfoListConstIter iter = m_typeInfo.begin(); iter != m_typeInfo.end(); iter++) {
    if (strType == (*iter).name) {
      return RbtTriposAtomType::eType(iter - m_typeInfo.begin());
    }
  }
  return RbtTriposAtomType::UNDEFINED;
}

//Setup atom type info
void RbtTriposAtomType::SetupTypeInfo() {
  m_typeInfo = RbtTriposInfoList(MAXTYPES,info());
  m_typeInfo[UNDEFINED] = info("UNDEFINED",0,RbtAtom::UNDEFINED);
  m_typeInfo[Al]        = info("Al",13,RbtAtom::UNDEFINED);
  m_typeInfo[Br]        = info("Br",35,RbtAtom::SP3);
  m_typeInfo[C_cat]     = info("C.cat",6,RbtAtom::SP2);
  m_typeInfo[C_1]       = info("C.1",6,RbtAtom::SP);
  m_typeInfo[C_1_H1]    = info("C.1.H1",6,RbtAtom::SP);
  m_typeInfo[C_2]       = info("C.2",6,RbtAtom::SP2);
  m_typeInfo[C_2_H1]    = info("C.2.H1",6,RbtAtom::SP2);
  m_typeInfo[C_2_H2]    = info("C.2.H2",6,RbtAtom::SP2);
  m_typeInfo[C_3]       = info("C.3",6,RbtAtom::SP3);
  m_typeInfo[C_3_H1]    = info("C.3.H1",6,RbtAtom::SP3);
  m_typeInfo[C_3_H2]    = info("C.3.H2",6,RbtAtom::SP3);
  m_typeInfo[C_3_H3]    = info("C.3.H3",6,RbtAtom::SP3);
  m_typeInfo[C_ar]      = info("C.ar",6,RbtAtom::AROM);
  m_typeInfo[C_ar_H1]   = info("C.ar.H1",6,RbtAtom::AROM);
  m_typeInfo[Ca]        = info("Ca",20,RbtAtom::UNDEFINED);
  m_typeInfo[Cl]        = info("Cl",17,RbtAtom::SP3);
  m_typeInfo[Du]        = info("Du",0,RbtAtom::UNDEFINED);
  m_typeInfo[F]         = info("F",9,RbtAtom::SP3);
  m_typeInfo[H]         = info("H",1,RbtAtom::SP3);
  m_typeInfo[H_P]       = info("H.P",1,RbtAtom::SP3);
  m_typeInfo[I]         = info("I",53,RbtAtom::SP3);
  m_typeInfo[K]         = info("K",19,RbtAtom::UNDEFINED);
  m_typeInfo[Li]        = info("Li",3,RbtAtom::UNDEFINED);
  m_typeInfo[LP]        = info("LP",0,RbtAtom::UNDEFINED);
  m_typeInfo[N_1]       = info("N.1",7,RbtAtom::SP);
  m_typeInfo[N_2]       = info("N.2",7,RbtAtom::SP2);
  m_typeInfo[N_3]       = info("N.3",7,RbtAtom::SP3);
  m_typeInfo[N_4]       = info("N.4",7,RbtAtom::SP3);
  m_typeInfo[N_am]      = info("N.am",7,RbtAtom::TRI);
  m_typeInfo[N_ar]      = info("N.ar",7,RbtAtom::AROM);
  m_typeInfo[N_pl3]     = info("N.pl3",7,RbtAtom::TRI);
  m_typeInfo[Na]        = info("Na",11,RbtAtom::UNDEFINED);
  m_typeInfo[O_co2]     = info("O.co2",8,RbtAtom::TRI);
  m_typeInfo[O_2]       = info("O.2",8,RbtAtom::SP2);
  m_typeInfo[O_3]       = info("O.3",8,RbtAtom::SP3);
  m_typeInfo[P_3]       = info("P.3",15,RbtAtom::SP3);
  m_typeInfo[S_2]       = info("S.2",16,RbtAtom::SP2);
  m_typeInfo[S_3]       = info("S.3",16,RbtAtom::SP3);
  m_typeInfo[S_o]       = info("S.o",16,RbtAtom::SP2);
  m_typeInfo[S_o2]      = info("S.o2",16,RbtAtom::SP2);
  m_typeInfo[Si]        = info("Si",14,RbtAtom::UNDEFINED);
}
