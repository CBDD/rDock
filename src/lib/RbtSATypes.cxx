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

#include "RbtSATypes.h"
#include "RbtAtom.h"
#include "RbtBond.h"

RbtHHSType::RbtHHSType()
{
  SetupTypeNames();
}

RbtHHSType::~RbtHHSType() {}


RbtHHSType::eType RbtHHSType::operator() (RbtAtom* anAtom) const {
  RbtAtomList bondedAtomList;
  RbtAtomPtr bondedAtom;
  RbtBool polar(false);
  Rbt::isAtomCationic bIsCat;
  Rbt::isAtomAnionic bIsAni;
  Rbt::isAtomMetal bIsMet;
  Rbt::isAtomLipophilic bIsLipo;
  RbtInt nH = anAtom->GetNumImplicitHydrogens() + anAtom->GetCoordinationNumber(1);
  switch (anAtom->GetAtomicNo()) {

  case 1://Hydrogen
    bondedAtomList = Rbt::GetBondedAtomList(anAtom);
    if (bondedAtomList.empty())
      return H;//Isolated hydrogen atom
    bondedAtom = bondedAtomList.front();
    switch (bondedAtom->GetAtomicNo()) {
    case 6: return H;
    case 7: return (bIsCat(anAtom)) ? HNp : HN;
    case 8: return HO;
    case 16: return HS;
    default: return UNDEFINED;
    }

  case 6://Carbon
    //Deal with special case of charged carbon (guanidinium, imidazolium) first
    //Declare as undefined if not SP2 (or AROM : DM 7 Jan 2004)
    if (bIsCat(anAtom)) {
      switch (anAtom->GetHybridState()) {
      case RbtAtom::SP2:
      case RbtAtom::AROM:
	return C_sp2p;
      default:
	return UNDEFINED;
      }
    }
    //Neutral carbons from here on
    //Loose definition of "polar" carbons: bonded to at least one O or N
    polar = ( (anAtom->GetCoordinationNumber(7) + anAtom->GetCoordinationNumber(8)) > 0);
    switch (anAtom->GetHybridState()) {
    case RbtAtom::SP3:
      switch (nH) {
      case 0:
	return (polar) ? C_sp3_P : C_sp3;
      case 1:
	return (polar) ? CH_sp3_P : CH_sp3;
      case 2:
	return (polar) ? CH2_sp3_P : CH2_sp3;
	//We don't have a type for CH4 (isolated methane) so double up with the CH3- params
	//Hardly likely to be used!
      case 3:
      case 4:
	return (polar) ? CH3_sp3_P : CH3_sp3;
      default:
	return UNDEFINED;
      }
    case RbtAtom::SP2:
      switch (nH) {
      case 0:
	return (polar) ? C_sp2_P : C_sp2;
      case 1:
	return (polar) ? CH_sp2_P : CH_sp2;
      case 2:
	return (polar) ? CH2_sp2_P : CH2_sp2;
      default:
	return UNDEFINED;
      }
    case RbtAtom::AROM:
      switch (nH) {
      case 0:
	return (polar) ? C_ar_P : C_ar;
      case 1:
	return (polar) ? CH_ar_P : CH_ar;
      default:
	return UNDEFINED;
      }
    case RbtAtom::SP:
      return C_sp;
    default:
      return UNDEFINED;
    }

  case 7://Nitrogen
    bondedAtomList = Rbt::GetBondedAtomList(anAtom);
    //Deal with charged nitrogens first. Defined as:
    //i) N carrying +ve charge, OR
    //ii) N bonded to 4 atoms, OR
    //iii) N bonded to +ve charged atom (more likely)
    if ( bIsCat(anAtom) || (bondedAtomList.size() == 4) ||
	 (Rbt::GetNumAtoms(bondedAtomList,bIsCat) > 0) ) {
      switch (anAtom->GetHybridState()) {
      case RbtAtom::SP3:
	return N_sp3p;
	//One of the charged nitrogens in a guanidinium/imidazolium will be
	//formally TRI: trigonal planar (two single bonds),
	//the other will be formally SP2 (one double/one single bond)
	//both should be defined as N_sp2p for the solvation calculation
      case RbtAtom::TRI:
      case RbtAtom::SP2:
      case RbtAtom::AROM:
	return N_sp2p;
      default:
	return UNDEFINED;
      }
    }
    //Neutral nitrogens from here on
    switch (anAtom->GetHybridState()) {
    case RbtAtom::SP3:
      switch (nH) {
      case 0:
	return N_sp3;
      case 1:
	return NH_sp3;
      case 2:
	return NH2_sp3;
      default:
	return UNDEFINED;
      }
    case RbtAtom::TRI:
      switch (nH) {
      case 0:
	return N_tri;
      case 1:
	return NH_tri;
      case 2:
	return NH2_tri;
      default:
	return UNDEFINED;
      }
    case RbtAtom::SP2:
      //N in nitro groups are formally defined as N_SP2+ by MdlFileSource
      //with 3 bonds, but have their formal charge neutralised by RemoveZwitterions
      //They need to be trapped here. Better to reuse the N_tri solvation type.
      return (anAtom->GetNumBonds() == 3) ? N_tri : N_sp2;
    case RbtAtom::AROM:
      return N_ar;
    case RbtAtom::SP:
      return N_sp;
    default:
      return UNDEFINED;
    }

  case 8: //Oxygen
    //Deal with negatively charged oxygens first
    if (bIsAni(anAtom))
      return Om;
    //Neutral oxygens from here on
    //Terminal oxygen bonded to nitrogen (mainly nitro groups)
    if ( (anAtom->GetCoordinationNumber(7) == 1) && (anAtom->GetNumBonds() == 1))
      return ON;
    switch (anAtom->GetHybridState()) {
    case RbtAtom::SP3:
      return (nH == 1) ? OH_sp3 : (nH == 2) ? OW : O_sp3;//Covers alcohols, ethers and isolated H2O
    case RbtAtom::TRI:
      //trap aromatic ethers and phenols here
      //MdlFileSource converts the O_SP3 to O_TRI but for solvation purposes
      //they are best considered as O_sp3
      //We check for O_TRI bonded to only lipophilic atoms or hydrogen
      bondedAtomList = Rbt::GetBondedAtomList(anAtom);
      if ( (Rbt::GetNumAtoms(bondedAtomList,bIsLipo) + nH) == 2)
	return (nH > 0) ? OH_sp3 : O_sp3;
      //O_tri is reserved for esters, OH_tri for acids
      return (nH > 0) ? OH_tri : O_tri;
    case RbtAtom::SP2:
      return O_sp2;
    default:
      return UNDEFINED;
    }

  case 16: //Sulphur
    switch (anAtom->GetHybridState()) {
    case RbtAtom::SP3:
    case RbtAtom::TRI:
      return S_sp3;
    case RbtAtom::SP2:
      return S_sp2;
    default:
      return UNDEFINED;
    }

  case 15: //Phosphorous
    return P;
  case 9: //Fluorine
    return F;
  case 17: //Chlorine
    return Cl;
  case 35: //Bromine
    return Br;
  case 53: //Iodine
    return I;

  default://Catch metals here, all else is undefined type
    return (bIsMet(anAtom)) ? Metal : UNDEFINED;
  }
}

RbtString RbtHHSType::Type2Str(RbtHHSType::eType aType) const {
  return m_typeNames[aType];
}

RbtHHSType::eType RbtHHSType::Str2Type(const RbtString& strType) const {
  RbtStringListConstIter iter = std::find(m_typeNames.begin(),m_typeNames.end(), strType);
  if (iter != m_typeNames.end()) {
    return RbtHHSType::eType(iter - m_typeNames.begin());
  }
  else {
    return RbtHHSType::UNDEFINED;
  }
}

void RbtHHSType::SetupTypeNames() {
  //Create a vector<string> of the right size
  m_typeNames = RbtStringList(RbtHHSType::MAXTYPES,RbtString());
  m_typeNames[RbtHHSType::UNDEFINED] = "UNDEFINED";
  m_typeNames[RbtHHSType::C_sp3]     = "C_sp3";
  m_typeNames[RbtHHSType::CH_sp3]    = "CH_sp3";
  m_typeNames[RbtHHSType::CH2_sp3]   = "CH2_sp3";
  m_typeNames[RbtHHSType::CH3_sp3]   = "CH3_sp3";
  m_typeNames[RbtHHSType::C_sp2]     = "C_sp2";
  m_typeNames[RbtHHSType::CH_sp2]    = "CH_sp2";
  m_typeNames[RbtHHSType::CH2_sp2]   = "CH2_sp2";
  m_typeNames[RbtHHSType::C_sp2p]    = "C_sp2p";
  m_typeNames[RbtHHSType::C_ar]      = "C_ar";
  m_typeNames[RbtHHSType::CH_ar]     = "CH_ar";
  m_typeNames[RbtHHSType::C_sp]      = "C_sp";
  m_typeNames[RbtHHSType::C_sp3_P]   = "C_sp3_P";
  m_typeNames[RbtHHSType::CH_sp3_P]  = "CH_sp3_P";
  m_typeNames[RbtHHSType::CH2_sp3_P] = "CH2_sp3_P";
  m_typeNames[RbtHHSType::CH3_sp3_P] = "CH3_sp3_P";
  m_typeNames[RbtHHSType::C_sp2_P]   = "C_sp2_P";
  m_typeNames[RbtHHSType::CH_sp2_P]  = "CH_sp2_P";
  m_typeNames[RbtHHSType::CH2_sp2_P] = "CH2_sp2_P";
  m_typeNames[RbtHHSType::C_ar_P]    = "C_ar_P";
  m_typeNames[RbtHHSType::CH_ar_P]   = "CH_ar_P";
  m_typeNames[RbtHHSType::H]         = "H";
  m_typeNames[RbtHHSType::HO]        = "HO";
  m_typeNames[RbtHHSType::HN]        = "HN";
  m_typeNames[RbtHHSType::HNp]       = "HNp";
  m_typeNames[RbtHHSType::HS]        = "HS";
  m_typeNames[RbtHHSType::O_sp3]     = "O_sp3";
  m_typeNames[RbtHHSType::OH_sp3]    = "OH_sp3";
  m_typeNames[RbtHHSType::OW]        = "OW";
  m_typeNames[RbtHHSType::O_tri]     = "O_tri";
  m_typeNames[RbtHHSType::OH_tri]    = "OH_tri";
  m_typeNames[RbtHHSType::O_sp2]     = "O_sp2";
  m_typeNames[RbtHHSType::ON]        = "ON";
  m_typeNames[RbtHHSType::Om]        = "Om";
  m_typeNames[RbtHHSType::N_sp3]     = "N_sp3";
  m_typeNames[RbtHHSType::NH_sp3]    = "NH_sp3";
  m_typeNames[RbtHHSType::NH2_sp3]   = "NH2_sp3";
  m_typeNames[RbtHHSType::N_sp3p]    = "N_sp3p";
  m_typeNames[RbtHHSType::N_tri]     = "N_tri";
  m_typeNames[RbtHHSType::NH_tri]    = "NH_tri";
  m_typeNames[RbtHHSType::NH2_tri]   = "NH2_tri";
  m_typeNames[RbtHHSType::N_sp2]     = "N_sp2";
  m_typeNames[RbtHHSType::N_sp2p]    = "N_sp2p";
  m_typeNames[RbtHHSType::N_ar]      = "N_ar";
  m_typeNames[RbtHHSType::N_sp]      = "N_sp";
  m_typeNames[RbtHHSType::S_sp3]     = "S_sp3";
  m_typeNames[RbtHHSType::S_sp2]     = "S_sp2";
  m_typeNames[RbtHHSType::P]         = "P";
  m_typeNames[RbtHHSType::F]         = "F";
  m_typeNames[RbtHHSType::Cl]        = "Cl";
  m_typeNames[RbtHHSType::Br]        = "Br";
  m_typeNames[RbtHHSType::I]         = "I";
  m_typeNames[RbtHHSType::Metal]     = "Metal";
}

//////////////////////
//Solvation class (extended atom)
//Constants from the HHS paper
const RbtDouble HHS_Solvation::r_s = 0.6;//solvent probe radius
const RbtDouble HHS_Solvation::d_s = 2.0 * HHS_Solvation::r_s;//solvent probe diameter
const RbtDouble HHS_Solvation::Pij_12 = 0.8875; // correction factor for 1-2 connections
const RbtDouble HHS_Solvation::Pij_13 = 0.3516; // ditto for 1-3
const RbtDouble HHS_Solvation::Pij_14 = 0.3156; // and for 1-4+

HHS_Solvation::HHS_Solvation(RbtHHSType::eType t, RbtAtom* a, RbtDouble p, RbtDouble r, RbtDouble s)
  : p_i(p), r_i(r), sigma(s), hhsType(t), atom(a) {
  Init();
}

//Actions
//Initialises A to 1.0, and S to area of isolated sphere
void HHS_Solvation::Init() {
  S_i = 4.0 * M_PI * (r_i + r_s) * (r_i + r_s);
  A_i = 1.0;
  A_inv = 1.0;
  E_i = S_i * sigma;
  m_var.clear();
  m_prt.clear();
  PI_r_i_plus_r_s = M_PI * (r_i + r_s);
  p_i_over_S_i = p_i / S_i;
}

//void HHS_Solvation::Save() {A_inv = A_i;}
//void HHS_Solvation::Restore() {A_i = A_inv;}

void HHS_Solvation::Overlap(HHS_Solvation* h, RbtDouble p_ij) {
  // check is there an overlap at all
  RbtDouble d2 = Rbt::Length2(atom->GetCoords(),h->atom->GetCoords());
  RbtDouble ol = r_i + h->r_i + 2.0*r_s;
  if( ol*ol < d2 ){
    return;
  }

  RbtDouble d = sqrt(d2);
  RbtDouble recip_d = 1.0 / d;
  RbtDouble ol_minus_d = ol - d;
  RbtDouble r_i_diff_over_d = recip_d * (h->r_i - r_i);
  //Update the area for this atom
  RbtDouble b_ij =  PI_r_i_plus_r_s * ol_minus_d * (1.0 + r_i_diff_over_d);
  RbtDouble A = 1.0 - (p_i_over_S_i * p_ij * b_ij);
  //Keep A between zero and one
  //otherwise get odd behavior for interpenetrating atoms
  //(possible at start of GA)
  if (A <= 0.0)
    A_i = 0.0;
  else if (A < 1.0)
    A_i *= A;
  //if A >= 1 do nothing (effectively A_i *= 1.0)
  
  //Update the area for the other atom (h)
  b_ij = h->PI_r_i_plus_r_s * ol_minus_d * (1.0 - r_i_diff_over_d);
  A = 1.0 - (h->p_i_over_S_i * p_ij * b_ij);
  if (A <= 0.0)
    h->A_i = 0.0;
  else if (A < 1.0)
    h->A_i *= A;
}

void HHS_Solvation::AddVariable(HHS_Solvation* anAtom) {
  m_var.push_back(anAtom);
  m_prt.push_back(anAtom);
}

//Helper function to calculate the overlap for all the variable interactions to this center
//All are 1-4+ by definition
void HHS_Solvation::OverlapVariable() {
  //Rbt::OverlapHHS overlap(this,HHS_Solvation::Pij_14);
  //std::for_each(m_prt.begin(),m_prt.end(),overlap);
  for (HHS_SolvationRListIter iter = m_prt.begin(); iter != m_prt.end(); ++iter)
    Overlap(*iter,HHS_Solvation::Pij_14);
}

//As above, but ignores disabled interaction centers
void HHS_Solvation::OverlapVariableEnabledOnly() {
  if (GetAtom()->GetEnabled()) {
    for (HHS_SolvationRListIter iter = m_prt.begin(); iter != m_prt.end(); ++iter) {
      if ((*iter)->GetAtom()->GetEnabled()) {
	Overlap(*iter,HHS_Solvation::Pij_14);
      }
    }
  }
}

//Partition the list of variable interactions
//Copy those that are closer than distance d to the partitioned list
//d=0 => clear the partition
void HHS_Solvation::Partition(RbtDouble d) {
  if (d > 0) {
    RbtDouble dd(d*d);
    m_prt.clear();
    for (HHS_SolvationRListConstIter iter = m_var.begin(); iter != m_var.end(); iter++) {
      if (Rbt::Length2((*iter)->atom->GetCoords(),atom->GetCoords()) < dd) {
	m_prt.push_back(*iter);
	//cout << "Partitioned " << (*iter)->atom->GetFullAtomName() << " - " << atom->GetFullAtomName() << endl;
      }
    }
  }
  else {
    m_prt = m_var;
  }
}

RbtBool Rbt::isHHSSelected::operator() (const HHS_Solvation* pHHS) const {
  return pHHS->GetAtom()->GetSelectionFlag();
}
