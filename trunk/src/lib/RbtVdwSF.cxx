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

#include "RbtVdwSF.h"
#include "RbtModel.h"

//Static data members
RbtString RbtVdwSF::_CT("RbtVdwSF");
RbtString RbtVdwSF::_USE_4_8("USE_4_8");
RbtString RbtVdwSF::_USE_TRIPOS("USE_TRIPOS");
RbtString RbtVdwSF::_RMAX("RMAX");
RbtString RbtVdwSF::_ECUT("ECUT");
RbtString RbtVdwSF::_E0("E0");

RbtVdwSF::RbtVdwSF() :
  m_use_4_8(true),m_use_tripos(false),m_rmax(1.5),m_ecut(1.0),m_e0(1.5)
{
#ifdef _DEBUG
  cout << _CT << " default constructor" << endl;
#endif //_DEBUG
  //Add parameters
  AddParameter(_USE_4_8,m_use_4_8);
  AddParameter(_USE_TRIPOS,m_use_tripos);
  AddParameter(_RMAX,m_rmax);
  AddParameter(_ECUT,m_ecut);
  AddParameter(_E0,m_e0);
  m_spVdwSource = RbtParameterFileSourcePtr(new RbtParameterFileSource(Rbt::GetRbtFileName("data/sf","Tripos52_vdw.prm")));
  Setup();
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtVdwSF::~RbtVdwSF() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}

//As this has a virtual base class we need a separate OwnParameterUpdated
//which can be called by concrete subclass ParameterUpdated methods
//See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
void RbtVdwSF::OwnParameterUpdated(const RbtString& strName) {
  //DM 25 Oct 2000 - heavily used params
  if (strName == _USE_4_8) {
    m_use_4_8 = GetParameter( _USE_4_8);
    Setup();
  }
  else if (strName == _USE_TRIPOS) {
    m_use_tripos = GetParameter(_USE_TRIPOS);
    Setup();
  }
  else if (strName == _RMAX) {
    m_rmax = GetParameter(_RMAX);
    Setup();
  }
  else if (strName == _ECUT) {
    m_ecut = GetParameter(_ECUT);
    SetupCloseRange();
  }
  else if (strName == _E0) {
    m_e0 = GetParameter(_E0);
    SetupCloseRange();
  }
}

//Used by subclasses to calculate vdW potential between pAtom and all atoms in atomList
RbtDouble RbtVdwSF::VdwScore(const RbtAtom* pAtom, const RbtAtomRList& atomList) const {
  RbtDouble score = 0.0;
  if (atomList.empty()) {
    return score;
  }

  const RbtCoord& c1 = pAtom->GetCoords();
  //Get the iterator into the appropriate row of the vdw table for this atom type
  RbtTriposAtomType::eType type1 = pAtom->GetTriposType();
  RbtVdwTableConstIter iter1 = m_vdwTable.begin() + type1;

  //4-8 potential, never annotated
  if (m_use_4_8) {
    for (RbtAtomRListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
      const RbtCoord& c2 = (*iter)->GetCoords();
      RbtDouble R_sq = Rbt::Length2(c1,c2);//Distance squared
      RbtTriposAtomType::eType type2 = (*iter)->GetTriposType();
      //iter2 points to the vdw params for this atom type pair
      RbtVdwRowConstIter iter2 = (*iter1).begin() + type2;
      RbtDouble s = f4_8(R_sq,*iter2);
      score += s;
    }
  }
  //6-12 with annotation
  else if (isAnnotationEnabled()) {
    for (RbtAtomRListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
      const RbtCoord& c2 = (*iter)->GetCoords();
      RbtDouble R_sq = Rbt::Length2(c1,c2);//Distance squared
      RbtTriposAtomType::eType type2 = (*iter)->GetTriposType();
      //iter2 points to the vdw params for this atom type pair
      RbtVdwRowConstIter iter2 = (*iter1).begin() + type2;
      RbtDouble s = f6_12(R_sq,*iter2);
      if (s != 0.0) {
	score += s;
	RbtAnnotationPtr spAnnotation(new RbtAnnotation(pAtom,*iter,sqrt(R_sq),s));
	AddAnnotation(spAnnotation);
      }
    }
  }
  //6-12 without annotation
  else {
    for (RbtAtomRListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
      const RbtCoord& c2 = (*iter)->GetCoords();
      RbtDouble R_sq = Rbt::Length2(c1,c2);//Distance squared
      RbtTriposAtomType::eType type2 = (*iter)->GetTriposType();
      //iter2 points to the vdw params for this atom type pair
      RbtVdwRowConstIter iter2 = (*iter1).begin() + type2;
      RbtDouble s = f6_12(R_sq,*iter2);
      score += s;
    }
  }
  return score;
}

//As above, but score is calculated only between enabled atoms
RbtDouble RbtVdwSF::VdwScoreEnabledOnly(const RbtAtom* pAtom, const RbtAtomRList& atomList) const {
  RbtDouble score = 0.0;
  if (!pAtom->GetEnabled() || atomList.empty()) {
    return score;
  }

  const RbtCoord& c1 = pAtom->GetCoords();
  //Get the iterator into the appropriate row of the vdw table for this atom type
  RbtTriposAtomType::eType type1 = pAtom->GetTriposType();
  RbtVdwTableConstIter iter1 = m_vdwTable.begin() + type1;

  //4-8 potential, never annotated
  if (m_use_4_8) {
    for (RbtAtomRListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
      if ((*iter)->GetEnabled()) {
	const RbtCoord& c2 = (*iter)->GetCoords();
	RbtDouble R_sq = Rbt::Length2(c1,c2);//Distance squared
	RbtTriposAtomType::eType type2 = (*iter)->GetTriposType();
	//iter2 points to the vdw params for this atom type pair
	RbtVdwRowConstIter iter2 = (*iter1).begin() + type2;
	RbtDouble s = f4_8(R_sq,*iter2);
	score += s;
      }
    }
  }
  //6-12 with annotation
  else if (isAnnotationEnabled()) {
    for (RbtAtomRListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
      if ((*iter)->GetEnabled()) {
	const RbtCoord& c2 = (*iter)->GetCoords();
	RbtDouble R_sq = Rbt::Length2(c1,c2);//Distance squared
	RbtTriposAtomType::eType type2 = (*iter)->GetTriposType();
	//iter2 points to the vdw params for this atom type pair
	RbtVdwRowConstIter iter2 = (*iter1).begin() + type2;
	RbtDouble s = f6_12(R_sq,*iter2);
	if (s != 0.0) {
	  score += s;
	  RbtAnnotationPtr spAnnotation(new RbtAnnotation(pAtom,*iter,sqrt(R_sq),s));
	  AddAnnotation(spAnnotation);
	}
      }
    }
  }
  //6-12 without annotation
  else {
    for (RbtAtomRListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
      if ((*iter)->GetEnabled()) {
	const RbtCoord& c2 = (*iter)->GetCoords();
	RbtDouble R_sq = Rbt::Length2(c1,c2);//Distance squared
	RbtTriposAtomType::eType type2 = (*iter)->GetTriposType();
	//iter2 points to the vdw params for this atom type pair
	RbtVdwRowConstIter iter2 = (*iter1).begin() + type2;
	RbtDouble s = f6_12(R_sq,*iter2);
	score += s;
      }
    }
  }
  return score;
}

RbtDouble RbtVdwSF::MaxVdwRange(const RbtAtom* pAtom) const {
  return m_maxRange[pAtom->GetTriposType()];
}

RbtDouble RbtVdwSF::MaxVdwRange(RbtTriposAtomType::eType t) const {
  return m_maxRange[t];
}

//Initialise m_vdwTable with appropriate params for each atom type pair
void RbtVdwSF::Setup() {
  RbtTriposAtomType triposType;
  RbtInt iTrace = GetTrace();
  if (iTrace > 3) {
    cout << endl << _CT << "::Setup()" << endl;
    cout << "TYPE1,TYPE2,A,B,RMIN,KIJ,RMAX" << endl;
  }
  //Dummy read to force parsing of file, otherwise the first SetSection is overridden
  RbtStringList secList = m_spVdwSource->GetSectionList();
  RbtString _R("R");
  RbtString _K("K");
  RbtString _IP("IP");
  RbtString _POL("POL");
  RbtString _ISHBD("isHBD");
  RbtString _ISHBA("isHBA");
  m_vdwTable = RbtVdwTable(RbtTriposAtomType::MAXTYPES,RbtVdwRow(RbtTriposAtomType::MAXTYPES));
  m_maxRange = RbtDoubleList(RbtTriposAtomType::MAXTYPES,0.0);
  for (RbtInt i = RbtTriposAtomType::UNDEFINED; i < RbtTriposAtomType::MAXTYPES; i++) {
    //Read the params for atom type i
    RbtString stri = triposType.Type2Str(RbtTriposAtomType::eType(i));
    m_spVdwSource->SetSection(stri);
    RbtDouble Ri = m_spVdwSource->GetParameterValue(_R);//vdw radius
    RbtDouble Ki = m_spVdwSource->GetParameterValue(_K);//Tripos 5.2 well depth
    RbtBool hasIPi = m_spVdwSource->isParameterPresent(_IP);
    RbtBool isHBDi = m_spVdwSource->isParameterPresent(_ISHBD);
    RbtBool isHBAi = m_spVdwSource->isParameterPresent(_ISHBA);
    RbtDouble Ii, alphai;
    if (hasIPi) {
      Ii = m_spVdwSource->GetParameterValue(_IP);//Ionisation potential, I
      alphai = m_spVdwSource->GetParameterValue(_POL);//Polarisability, alpha
    }
    //Symmetric matrix, so only need to generate one half (i.e. m_vdwTable[i][j] = m_vdwTable[j][i])
    for (RbtInt j = i; j < RbtTriposAtomType::MAXTYPES; j++) {
      //Read the params for atom type j
      RbtString strj = triposType.Type2Str(RbtTriposAtomType::eType(j));
      m_spVdwSource->SetSection(strj);
      RbtDouble Rj = m_spVdwSource->GetParameterValue(_R);//vdw radius
      RbtDouble Kj = m_spVdwSource->GetParameterValue(_K);//Tripos 5.2 well depth
      RbtBool hasIPj = m_spVdwSource->isParameterPresent(_IP);
      RbtBool isHBDj = m_spVdwSource->isParameterPresent(_ISHBD);
      RbtBool isHBAj = m_spVdwSource->isParameterPresent(_ISHBA);
      RbtDouble Ij, alphaj;
      if (hasIPj) {
	Ij = m_spVdwSource->GetParameterValue(_IP);//Ionisation potential, I
	alphaj = m_spVdwSource->GetParameterValue(_POL);//Polarisability, alpha
      }

      //Create the vdw params for this atom type pair
      vdwprms prms;
      prms.rmin = Ri+Rj;//rmin is the sum of Ri and Rj
      RbtDouble rmax = prms.rmin * m_rmax;
      m_maxRange[i] = std::max(m_maxRange[i],rmax);//Keep track of max range for atom type i
      m_maxRange[j] = std::max(m_maxRange[j],rmax);//Keep track of max range for atom type j
      prms.rmax_sq = rmax * rmax;//Max range**2
      //EITHER: Well depth is zero between donor Hs and acceptors
      if ( (isHBDi && isHBAj) || (isHBDj && isHBAi)) {
	prms.kij = 0.0;
      }
      //OR, use Standard L-J combination rules for well depth
      //switch to this mode if either atom is missing IP values
      else if (m_use_tripos || !hasIPi || !hasIPj) {
	prms.kij = sqrt(Ki*Kj);
      }
      //OR, use GOLD rules based on ionisation potential and polarisability
      else {
	RbtDouble D = 0.345 * Ii * Ij * alphai * alphaj / (Ii + Ij);
	RbtDouble C = 0.5 * D * pow(prms.rmin,6);
	prms.kij = D * D / (4.0 * C);
      }
      //Having got the well depth, we can now generate A and B for either 4-8 or 6-12
      RbtDouble rmin_pwr = (m_use_4_8) ? pow(prms.rmin,4) : pow(prms.rmin,6);
      prms.A = prms.kij * rmin_pwr * rmin_pwr;
      prms.B = 2.0 * prms.kij * rmin_pwr;
      m_vdwTable[i][j] = prms;
      m_vdwTable[j][i] = prms;
      if (iTrace > 3) {
	cout << triposType.Type2Str(RbtTriposAtomType::eType(i)) << ","
	     << triposType.Type2Str(RbtTriposAtomType::eType(j)) << ","
	     << prms.A << ","
	     << prms.B << ","
	     << prms.rmin << ","
	     << prms.kij << ","
	     << sqrt(prms.rmax_sq) << endl;
      }
    }
  }
  //Now we can regenerate the close range params
  SetupCloseRange();
  //Set the overall range automatically from the max range for each atom type
  RbtDouble range = *(std::max_element(m_maxRange.begin(),m_maxRange.end()));
  SetRange(range);
  if (iTrace > 1) {
    cout << "Overall max range for scoring function = " << range << endl;
  }
}

//Regenerate the short-range params only (called more frequently)
void RbtVdwSF::SetupCloseRange() {
  RbtTriposAtomType triposType;
  RbtInt iTrace = GetTrace();
  if (iTrace > 3) {
    cout << endl << _CT << "::SetupCloseRange()" << endl;
    cout << "TYPE1,TYPE2,RCUT,ECUT,SLOPE,E0" << endl;
  }
  RbtDouble x = 1.0 + sqrt(1.0 + m_ecut);
  RbtDouble p = (m_use_4_8) ? pow(x,1.0/4.0) : pow(x,1.0/6.0);
  RbtDouble c = 1.0/p;
  for (RbtVdwTableIter iter1 = m_vdwTable.begin(); iter1 != m_vdwTable.end(); iter1++) {
    for (RbtVdwRowIter iter2 = (*iter1).begin(); iter2 != (*iter1).end(); iter2++) {
      (*iter2).rcutoff_sq = pow((*iter2).rmin * c, 2);
      (*iter2).ecutoff = (*iter2).kij * m_ecut;
      (*iter2).e0 = (*iter2).ecutoff * m_e0;
      (*iter2).slope = ((*iter2).e0 - (*iter2).ecutoff) / (*iter2).rcutoff_sq;
      if (iTrace > 3) {
	cout << triposType.Type2Str(RbtTriposAtomType::eType(iter1-m_vdwTable.begin())) << ","
	     << triposType.Type2Str(RbtTriposAtomType::eType(iter2-(*iter1).begin())) << ","
	     << sqrt((*iter2).rcutoff_sq) << ","
	     << (*iter2).ecutoff << ","
	     << (*iter2).slope << ","
	     << (*iter2).e0 << endl;
      }
    }
  }
}

//Index the flexible interactions between two atom lists.
//Foreach atom in the first list, compile a list of the atoms in the second list whose distance can vary to that atom
//intns container should have been initialised to a vector of empty atom lists prior to calling BuildIntraMap
void RbtVdwSF::BuildIntraMap(const RbtAtomRList& atomList1, const RbtAtomRList& atomList2, RbtAtomRListList& intns) const {
  Rbt::SelectAtom selectAtom(true);
  Rbt::isAtomSelected isSelected;
  RbtBool bSingleList = atomList2.empty();//If true, then index the flexible interactions within atomList1
  for (RbtAtomRListConstIter iter = atomList1.begin(); iter != atomList1.end(); iter++) {
  	//1. First select all atoms in the list
  	std::for_each(atomList1.begin(),atomList1.end(),selectAtom);
    RbtAtom* pAtom = (*iter);
    RbtInt id = pAtom->GetAtomId()-1;
    Assert<RbtAssert>(id >= 0 && id < intns.size());
    RbtModel* pModel = pAtom->GetModelPtr();
    //2. Now deselect all atoms in the same model as the current atom
    pModel->SetAtomSelectionFlags(false);
    //3. Now select all the flexible interactions to this atom within the current model
    //(Atoms in different models are not affected, therefore remain selected from step 1)
    //This way, we deal correctly with situations where the atoms in atomList may belong to different models
    pModel->SelectFlexAtoms(pAtom);
    if (bSingleList) {
    	std::copy_if(iter+1,atomList1.end(),std::back_inserter(intns[id]),isSelected);
    }
    else {
    	std::copy_if(atomList2.begin(),atomList2.end(),std::back_inserter(intns[id]),isSelected);
    }
    if (GetTrace() > 2) {
      cout << _CT << "::BuildIntraMap: " << pAtom->GetFullAtomName();
      if (!intns[id].empty()) {
	cout << "\t" << intns[id].front()->GetFullAtomName()
	     << " to " << intns[id].back()->GetFullAtomName();
      }
      cout << "; N = " << intns[id].size() << endl;
    }
  }
}

//Convenience method for above.
//Indexes the flexible interactions within a single atom list
void RbtVdwSF::BuildIntraMap(const RbtAtomRList& atomList, RbtAtomRListList& intns) const {
	RbtAtomRList emptyList;
	BuildIntraMap(atomList,emptyList,intns);
}


void RbtVdwSF::Partition(const RbtAtomRList& atomList, const RbtAtomRListList& intns, RbtAtomRListList& prtIntns, RbtDouble dist) const {
  RbtInt iTrace = GetTrace();

  //Clear the existing partitioned lists
  for (RbtAtomRListListIter iter = prtIntns.begin(); iter != prtIntns.end() ; iter++)
    (*iter).clear();
  if (iTrace > 3) {
    cout << _CT << ": Partition (dist=" << dist << ")" << endl;
  }

  for (RbtAtomRListConstIter iter = atomList.begin(); iter != atomList.end(); iter++) {
    RbtInt id = (*iter)->GetAtomId()-1;
    if (dist > 0.0) {
      isD_lt bInRange(*iter,dist);
      //This copies all interactions which are within dist A of atom
      std::copy_if(intns[id].begin(),intns[id].end(),std::back_inserter(prtIntns[id]),bInRange);
    }
    else {
      //Remove partition - simply copy from the master list of vdW interactions
      prtIntns[id] = intns[id];
    }
    if (iTrace > 3) {
      cout << _CT << ": " << (*iter)->GetFullAtomName() << ": #vdw=" << intns[id].size()
	   << ": #prtn=" << prtIntns[id].size() << endl;
    }
  }
}
