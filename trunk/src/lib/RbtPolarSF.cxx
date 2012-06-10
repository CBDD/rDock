/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtPolarSF.cxx#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtPolarSF.h"
#include "RbtWorkSpace.h"
#include "RbtPlane.h"

//Static data members
RbtString RbtPolarSF::_CT("RbtPolarSF");
RbtString RbtPolarSF::_R12FACTOR("R12FACTOR");
RbtString RbtPolarSF::_R12INCR("R12INCR");
RbtString RbtPolarSF::_DR12MIN("DR12MIN");
RbtString RbtPolarSF::_DR12MAX("DR12MAX");
RbtString RbtPolarSF::_A1("A1");
RbtString RbtPolarSF::_DA1MIN("DA1MIN");
RbtString RbtPolarSF::_DA1MAX("DA1MAX");
RbtString RbtPolarSF::_A2("A2");
RbtString RbtPolarSF::_DA2MIN("DA2MIN");
RbtString RbtPolarSF::_DA2MAX("DA2MAX");
RbtString RbtPolarSF::_INCMETAL("INCMETAL");
RbtString RbtPolarSF::_INCHBD("INCHBD");
RbtString RbtPolarSF::_INCHBA("INCHBA");
RbtString RbtPolarSF::_INCGUAN("INCGUAN");
RbtString RbtPolarSF::_GUAN_PLANE("GUAN_PLANE");
RbtString RbtPolarSF::_ABS_DR12("ABS_DR12");
RbtString RbtPolarSF::_LP_OSP2("LP_OSP2");
RbtString RbtPolarSF::_LP_PHI("LP_PHI");
RbtString RbtPolarSF::_LP_DPHIMIN("LP_DPHIMIN");
RbtString RbtPolarSF::_LP_DPHIMAX("LP_DPHIMAX");
RbtString RbtPolarSF::_LP_DTHETAMIN("LP_DTHETAMIN");
RbtString RbtPolarSF::_LP_DTHETAMAX("LP_DTHETAMAX");

RbtPolarSF::RbtPolarSF() :
  m_R12Factor(1.0),m_R12Incr(0.6),m_DR12Min(0.25),m_DR12Max(0.6),
  m_A1(180.0),m_DA1Min(30.0),m_DA1Max(80.0),m_A2(150.0),m_DA2Min(30.0),m_DA2Max(70.0),
  m_bAbsDR12(true),
  m_LP_PHI(45.0),m_LP_DPHIMin(15.0),m_LP_DPHIMax(30.0),
  m_LP_DTHETAMin(20.0),m_LP_DTHETAMax(60.0),
  m_PHI_lp_prms(m_LP_PHI,m_LP_DPHIMin,m_LP_DPHIMax),
  m_PHI_plane_prms(0.0,m_LP_PHI+m_LP_DPHIMin,m_LP_PHI+m_LP_DPHIMax),
  m_THETAprms(0.0,m_LP_DTHETAMin,m_LP_DTHETAMax)
{
#ifdef _DEBUG
  cout << _CT << " default constructor" << endl;
#endif //_DEBUG
  //Add parameters
  AddParameter(_R12FACTOR,m_R12Factor);
  AddParameter(_R12INCR,m_R12Incr);
  AddParameter(_DR12MIN,m_DR12Min);
  AddParameter(_DR12MAX,m_DR12Max);
  AddParameter(_A1,m_A1);
  AddParameter(_DA1MIN,m_DA1Min);
  AddParameter(_DA1MAX,m_DA1Max);
  AddParameter(_A2,m_A2);
  AddParameter(_DA2MIN,m_DA2Min);
  AddParameter(_DA2MAX,m_DA2Max);
  AddParameter(_INCMETAL,true);
  AddParameter(_INCHBD,true);
  AddParameter(_INCHBA,true);
  AddParameter(_INCGUAN,true);
  AddParameter(_GUAN_PLANE,true);
  AddParameter(_ABS_DR12,m_bAbsDR12);
  AddParameter(_LP_OSP2,false);
  AddParameter(_LP_PHI,m_LP_PHI);
  AddParameter(_LP_DPHIMIN,m_LP_DPHIMin);
  AddParameter(_LP_DPHIMAX,m_LP_DPHIMax);
  AddParameter(_LP_DTHETAMIN,m_LP_DTHETAMin);
  AddParameter(_LP_DTHETAMAX,m_LP_DTHETAMax);
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtPolarSF::~RbtPolarSF() {
#ifdef _DEBUG
  cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
	

RbtInteractionCenterList RbtPolarSF::CreateDonorInteractionCenters(const RbtAtomList& atomList) const
{
  const RbtBool bIncHBD = GetParameter(_INCHBD);
  const RbtBool bIncMetal = GetParameter(_INCMETAL);
  const RbtBool bIncGuan = GetParameter(_INCGUAN);
  const RbtBool bGuanPlane = GetParameter(_GUAN_PLANE);

  RbtInteractionCenterList intnList;

  //H-Bond donors
  if (bIncHBD) {
    RbtAtomList posList = Rbt::GetAtomList(atomList,Rbt::isAtomHBondDonor());
    for (RbtAtomListConstIter iter = posList.begin(); iter != posList.end(); iter++) {
      RbtAtomList parentList = Rbt::GetBondedAtomList(*iter);
      //Create interaction center from donor H and donor parent
      intnList.push_back(new RbtInteractionCenter(*iter,parentList.front()));
    }
  }

  //Metals
  if (bIncMetal) {
    RbtAtomList posList = Rbt::GetAtomList(atomList,Rbt::isAtomMetal());
    for (RbtAtomListConstIter iter = posList.begin(); iter != posList.end(); iter++) {
      //Create interaction center from metal ion only
      intnList.push_back(new RbtInteractionCenter(*iter));
    }
  }
  
  //Guanidinium carbons
  if (bIncGuan) {
    RbtAtomList posList = Rbt::GetAtomList(atomList,Rbt::isAtomGuanidiniumCarbon());
    for (RbtAtomListConstIter iter = posList.begin(); iter != posList.end(); iter++) {
      RbtAtomList parentList = Rbt::GetBondedAtomList(*iter);
      //Two options: for attractive interactions (HBA..GUAN) the SF depends on the angle between
      //the acceptor and the normal to the plane of the guanidinium, so we define an interaction
      //center with 3 atoms
      if (bGuanPlane) {
	//Create interaction center from guanidinium carbon and first 2 bonded atoms
	intnList.push_back(new RbtInteractionCenter(*iter,parentList[0],parentList[1]));
      }
      //For repulsive interactions (HBD..GUAN or GUAN..GUAN) the SF depends only on distance,
      //so we define an interaction center with just 1 atom
      else {
	//Create interaction center from guanidinium carbon only
	intnList.push_back(new RbtInteractionCenter(*iter));
      }	
    }
  }

  return intnList;
}

RbtInteractionCenterList RbtPolarSF::CreateAcceptorInteractionCenters(const RbtAtomList& atomList) const
{
  const RbtBool bIncHBA = GetParameter(_INCHBA);
  const RbtBool bLP = GetParameter(_LP_OSP2);
  
  RbtInteractionCenterList intnList;
  Rbt::isAtomicNo_eq bIsC(6);
  Rbt::isAtomicNo_eq bIsN(7);
  Rbt::isAtomicNo_eq bIsO(8);
  Rbt::isAtomRNA bIsRNA;
  Rbt::isAtomAnionic bIsAnionic;

  //Negative polar (HBA)
  if (bIncHBA) {
    RbtAtomList negList = Rbt::GetAtomList(atomList,Rbt::isAtomHBondAcceptor());
    for (RbtAtomListConstIter iter = negList.begin(); iter != negList.end(); iter++) {
      RbtAtomList parentList = Rbt::GetBondedAtomList(*iter);

      //Unconnected acceptor - ignore
      if (parentList.empty()) continue;

      //Singly connected acceptor
      if (parentList.size() == 1) {
	RbtAtomPtr spAcceptorParent = parentList.front();
	//Detect terminal oxygen bonded to carbon or nitrogen here
	//and define an interaction center including lone pair geometry if required
	if (bLP && bIsO(*iter) && (bIsC(spAcceptorParent) || bIsN(spAcceptorParent))) {
	  //Get the atoms bonded to the acceptor parent (not including the acceptor itself)
	  RbtAtomList grandParentList = Rbt::GetAtomList(Rbt::GetBondedAtomList(spAcceptorParent),
							 std::not1(std::bind2nd(Rbt::isAtomPtr_eq(),*iter)));
	  if (!grandParentList.empty()) {
	    RbtAtomPtr spGrandParent = grandParentList.front();
	    //If O is -ve charged (likely carboxylate), or is in RNA, or is in nitro group
	    //define interaction center with LONEPAIR geometry (directionality and in-plane of LP)
	    if (bIsAnionic(*iter) || bIsRNA(*iter) || bIsN(spAcceptorParent)) {
	      intnList.push_back(new RbtInteractionCenter(*iter,spAcceptorParent,spGrandParent,RbtInteractionCenter::LONEPAIR));
	      //cout << "LONEPAIR acceptor: " << (*iter)->GetFullAtomName() << endl;
	    }
	    //else define interaction center with PLANE geometry (in-plane of LP, broader directionality)
	    else {
	      intnList.push_back(new RbtInteractionCenter(*iter,spAcceptorParent,spGrandParent,RbtInteractionCenter::PLANE));
	      //cout << "PLANE acceptor: " << (*iter)->GetFullAtomName() << endl;
	    }	      
	  }
	  else {
	    intnList.push_back(new RbtInteractionCenter(*iter,spAcceptorParent));
	  }
	}
	else {
	  intnList.push_back(new RbtInteractionCenter(*iter,spAcceptorParent));
	}
      }
      else {
	//define pseudo atom for acceptor parent if more than one parent
    	//pseudo atoms stored permanently in receptor model
	RbtPseudoAtomPtr spPseudoAtom = (*iter)->GetModelPtr()->AddPseudoAtom(parentList);
	RbtAtomPtr spAcceptorParent = spPseudoAtom;
	intnList.push_back(new RbtInteractionCenter(*iter,spAcceptorParent));
      }
    }
  }
  return intnList;
}

void RbtPolarSF::BuildIntraMap(const RbtInteractionCenterList& ICList1,
			       const RbtInteractionCenterList& ICList2,
			       RbtInteractionListMap& intns) const {
  Rbt::SelectInteractionCenter selectIC(true);
  Rbt::isInteractionCenterSelected isSelected;
  RbtBool bSingleList = ICList2.empty();//If true, then index the flexible interactions within ICList1
  for (RbtInteractionCenterListConstIter iter = ICList1.begin(); iter != ICList1.end(); iter++) {
  	//1. First select all atoms of all the interaction centers in the list
  	//(Needed to deal with situations where the interaction centers do not all belong to the same model)
  	std::for_each(ICList1.begin(),ICList1.end(),selectIC);
    //We use the ID of Atom1 of the interaction center (IC) as the index into the map
    //but we check all the constituent atoms of the IC to assess whether the interaction should be deemed "flexible"
    //Reason is that some interactions have an angular dependence: the primary distance may be rigid, but the angular
    //dependence may not be
    RbtAtom* pAtom = (*iter)->GetAtom1Ptr();
    RbtInt id = pAtom->GetAtomId()-1;
    Assert<RbtAssert>(id >= 0 && id < intns.size());
    //2. Deselect all atoms in the same model as the current atom
    pAtom->GetModelPtr()->SetAtomSelectionFlags(false);
    //3. Set the selection flags for all atoms in the model whose distance can vary to *any* of the
    //atoms in the IC
    RbtAtomRList atomList = (*iter)->GetAtomList();
    std::for_each(atomList.begin(),atomList.end(),Rbt::SelectFlexAtoms());
    //Now we retrieve the ICs from the 2nd list whose constituent atoms have been selected
    if (bSingleList) {
    	std::copy_if(iter+1,ICList1.end(),std::back_inserter(intns[id]),isSelected);
    }
    else {
    	std::copy_if(ICList2.begin(),ICList2.end(),std::back_inserter(intns[id]),isSelected);
    }
    if (GetTrace() > 2) {
      cout << GetFullName() << ": Intns to " << pAtom->GetFullAtomName() << endl;
      for (RbtInteractionCenterListConstIter i = intns[id].begin(); i != intns[id].end(); i++) {
		cout << "\t" << (*i)->GetAtom1Ptr()->GetFullAtomName() << endl;
      }
    }
  }
}

//Convenience method for above.
//Indexes the flexible interactions within a single interaction center list.
void RbtPolarSF::BuildIntraMap(const RbtInteractionCenterList& ICList,
			       RbtInteractionListMap& intns) const {
	RbtInteractionCenterList emptyList;
	BuildIntraMap(ICList,emptyList,intns);
}

RbtDouble RbtPolarSF::IntraScore(const RbtInteractionCenterList& posList,
				 const RbtInteractionCenterList& negList,
				 const RbtInteractionListMap& intns, RbtBool attr) const {
  RbtDouble score = 0.0;//Total score

  RbtPolarSF::f1prms Rprms = GetRprms();//Distance params
  RbtPolarSF::f1prms A1prms = GetA1prms();//Donor angle params
  RbtPolarSF::f1prms A2prms = GetA2prms();//Acceptor angle params

  RbtDouble s(0.0);//Partial scores

  if (attr) {
    //Pos-neg
    for (RbtInteractionCenterListConstIter iter = posList.begin(); iter != posList.end(); iter++) {
      RbtAtom* pAtom = (*iter)->GetAtom1Ptr();
      RbtInt id = pAtom->GetAtomId()-1;
      //NO CHECK ON ID IN RANGE
      s = PolarScore(*iter,intns[id],Rprms,A1prms,A2prms);
      //if (fabs(s) > 0) {
      //cout << GetFullName() << ": pos-neg " << pAtom->GetFullAtomName() << "; s=" << s << endl;
      //}
      score += pAtom->GetUser1Value()*s;
    }
    //Neg-pos
    for (RbtInteractionCenterListConstIter iter = negList.begin(); iter != negList.end(); iter++) {
      RbtAtom* pAtom = (*iter)->GetAtom1Ptr();
      RbtInt id = pAtom->GetAtomId()-1;
      //NO CHECK ON ID IN RANGE
      s = PolarScore(*iter,intns[id],Rprms,A2prms,A1prms);
      //if (fabs(s) > 0) {
      //cout << GetFullName() << ": neg-pos " << pAtom->GetFullAtomName() << "; s=" << s << endl;
      //}
      score += pAtom->GetUser1Value()*s;
    }
  }
  else {
    //pos-pos
    for (RbtInteractionCenterListConstIter iter = posList.begin(); iter != posList.end(); iter++) {
      RbtAtom* pAtom = (*iter)->GetAtom1Ptr();
      RbtInt id = pAtom->GetAtomId()-1;
      //NO CHECK ON ID IN RANGE
      s = PolarScore(*iter,intns[id],Rprms,A1prms,A1prms);
      //if (fabs(s) > 0) {
      //	cout << GetFullName() << ": pos-pos " << pAtom->GetFullAtomName() << "; s=" << s << endl;
      //}
      score += pAtom->GetUser1Value()*s;
    }
    //neg-neg
    for (RbtInteractionCenterListConstIter iter = negList.begin(); iter != negList.end(); iter++) {
      RbtAtom* pAtom = (*iter)->GetAtom1Ptr();
      RbtInt id = pAtom->GetAtomId()-1;
      //NO CHECK ON ID IN RANGE
      s = PolarScore(*iter,intns[id],Rprms,A2prms,A2prms);
      //if (fabs(s) > 0) {
      //	cout << GetFullName() << ": neg-neg " << pAtom->GetFullAtomName() << "; s=" << s << endl;
      //}
      score += pAtom->GetUser1Value()*s;
    }
  }
  return score;
}

void RbtPolarSF::Partition(const RbtInteractionCenterList& posList,
			   const RbtInteractionCenterList& negList,
			   const RbtInteractionListMap& intns,
			   RbtInteractionListMap& prtIntns, RbtDouble dist) const {
  RbtInt iTrace = GetTrace();

  //Clear the existing partitioned lists
  for (RbtInteractionListMapIter iter = prtIntns.begin(); iter != prtIntns.end() ; iter++)
    (*iter).clear();
  if (iTrace > 3) {
    cout << "Partition (dist=" << dist << ")" << endl;
  }

  for (RbtInteractionCenterListConstIter iter = posList.begin(); iter != posList.end(); iter++) {
    RbtAtom* pAtom = (*iter)->GetAtom1Ptr();
    RbtInt id = pAtom->GetAtomId()-1;
    if (dist > 0.0) {
      Rbt::isInteractionD_lt bInRange(*iter,dist);
      //This copies all interactions which are within dist A of atom
      std::copy_if(intns[id].begin(),intns[id].end(),std::back_inserter(prtIntns[id]),bInRange);
    }
    else {
      //Remove partition - simply copy from the master list of vdW interactions
      prtIntns[id] = intns[id];
    }
    if (iTrace > 3) {
      cout << (*iter)->GetAtom1Ptr()->GetFullAtomName() << ": #intn=" << intns[id].size()
	   << ": #prtn=" << prtIntns[id].size() << endl;
    }
  }

  for (RbtInteractionCenterListConstIter iter = negList.begin(); iter != negList.end(); iter++) {
    RbtAtom* pAtom = (*iter)->GetAtom1Ptr();
    RbtInt id = pAtom->GetAtomId()-1;
    if (dist > 0.0) {
      Rbt::isInteractionD_lt bInRange(*iter,dist);
      //This copies all interactions which are within dist A of atom
      std::copy_if(intns[id].begin(),intns[id].end(),std::back_inserter(prtIntns[id]),bInRange);
    }
    else {
      //Remove partition - simply copy from the master list of vdW interactions
      prtIntns[id] = intns[id];
    }
    if (iTrace > 3) {
      cout << (*iter)->GetAtom1Ptr()->GetFullAtomName() << ": #intn=" << intns[id].size()
	   << ": #prtn=" << prtIntns[id].size() << endl;
    }
  }
}

RbtDouble RbtPolarSF::PolarScore(const RbtInteractionCenter* pIC1, const RbtInteractionCenterList& IC2List,
				 const f1prms& Rprms, const f1prms& A1prms, const f1prms& A2prms) const
{
  RbtDouble s(0.0);
  if (IC2List.empty()) {
    return s;
  }

  RbtCoord nullCoord;
  RbtBool bAnnotate = isAnnotationEnabled();
  //RbtString strName = GetName();

  RbtAtom* pAtom1_1 = pIC1->GetAtom1Ptr();
  //DM 7 June 2006 - check for disabled interaction centre 1
  if (!pAtom1_1->GetEnabled()) {
    return s;
  }
  const RbtCoord& cAtom1_1 = pAtom1_1->GetCoords();
  RbtAtom* pAtom1_2 = pIC1->GetAtom2Ptr();
  RbtAtom* pAtom1_3 = pIC1->GetAtom3Ptr();
  RbtInteractionCenter::eLP eLP1 = pIC1->LP();
  RbtBool bAngle1 = ( (pAtom1_2 != NULL) && (pAtom1_3 == NULL) ) ? true : false;
  RbtBool bPlane1 = ( (pAtom1_2 != NULL) && (pAtom1_3 != NULL) && (eLP1 == RbtInteractionCenter::NONE)) ? true : false;
  RbtBool bLP1 = ( (pAtom1_2 != NULL) && (pAtom1_3 != NULL) && (eLP1 != RbtInteractionCenter::NONE)) ? true : false;
  const f1prms& PHI1prms = (eLP1 == RbtInteractionCenter::LONEPAIR) ? m_PHI_lp_prms : m_PHI_plane_prms;
  const RbtCoord& cAtom1_2 = (bAngle1 || bPlane1 || bLP1) ? pAtom1_2->GetCoords(): nullCoord;
  const RbtCoord& cAtom1_3 = (bPlane1 || bLP1) ? pAtom1_3->GetCoords(): nullCoord;
  RbtPlane pl1 = (bPlane1 || bLP1) ? RbtPlane(cAtom1_1,cAtom1_2,cAtom1_3) : RbtPlane();
  RbtDouble radius1 = pAtom1_1->GetVdwRadius();

  for (RbtInteractionCenterListConstIter IC2Iter = IC2List.begin(); IC2Iter != IC2List.end(); IC2Iter++) {
    RbtAtom* pAtom2_1 = (*IC2Iter)->GetAtom1Ptr();
    //if (pAtom1_1 == pAtom2_1) continue;//check for self-interactions
    if (!pAtom2_1->GetEnabled()) continue;//check for disabled interaction centre 2
    const RbtCoord& cAtom2_1 = pAtom2_1->GetCoords();
    RbtAtom* pAtom2_2 = (*IC2Iter)->GetAtom2Ptr();
    RbtAtom* pAtom2_3 = (*IC2Iter)->GetAtom3Ptr();
    RbtInteractionCenter::eLP eLP2 = (*IC2Iter)->LP();
    RbtBool bAngle2 = ( (pAtom2_2 != NULL) && (pAtom2_3 == NULL) ) ? true : false;
    RbtBool bPlane2 = ( (pAtom2_2 != NULL) && (pAtom2_3 != NULL) && (eLP2 == RbtInteractionCenter::NONE)) ? true : false;
    RbtBool bLP2 = ( (pAtom2_2 != NULL) && (pAtom2_3 != NULL) && (eLP2 != RbtInteractionCenter::NONE)) ? true : false;
    RbtDouble radius2 = pAtom2_1->GetVdwRadius();
    RbtDouble R12 = m_R12Factor*(radius1+radius2)+m_R12Incr;
    RbtVector v12 = cAtom1_1 - cAtom2_1;
    RbtDouble R = v12.Length();
    RbtDouble DR = R-R12;
    RbtDouble f = m_bAbsDR12 ? f1(fabs(DR),Rprms) : f1(DR,Rprms);
    if (f > 0.0) {
      //cout << endl << "R12:\t" << pAtom1_1->GetFullAtomName() << "," << pAtom2_1->GetFullAtomName() << "," << R12
      //	   << "," << DR << ", f=" << f << endl;
      //For guanidinium bPlane2 interacting with C=O lone pair bLP1, we want to use the regular angular dependence
      if (bAngle1 || (bPlane2 && bLP1)) {
	RbtDouble DA1 = Rbt::Angle(cAtom1_2,cAtom1_1,cAtom2_1)-A1prms.R0;
	f *= f1(fabs(DA1),A1prms);
	//cout << "A1:\t" << A1prms.R0 << "," << DA1 << ", f=" << f << endl;
      }
      else if (bPlane1) {
	RbtDouble A = acos(-fabs(Rbt::Dot(v12.Unit(),pl1.VNorm())))*180.0/M_PI;
	RbtDouble DA1 = A-A1prms.R0;
	f *= f1(fabs(DA1),A1prms);
	//cout << "Pl1:\t" << A1prms.R0 << "," << DA1 << ", f=" << f << endl;
      }
      //Lone pair geometry dependence around Osp2
      else if (bLP1) {
	//Perpendicular distance from donor H to plane of lone pairs
	RbtDouble dPerp = Rbt::DistanceFromPointToPlane(cAtom2_1, pl1);
	//Coordinate of donor H projected into plane of lone pairs
	RbtCoord cPerp = cAtom2_1 - dPerp * pl1.VNorm();
	//cout << "dPerp = " << dPerp << " check cPerp = " << Rbt::DistanceFromPointToPlane(cPerp, pl1) << endl;
	//Can calculate sin(theta) directly from the two distances we have already
	RbtDouble theta = asin(dPerp/R)*180.0/M_PI;
	f *= f1(fabs(theta),m_THETAprms);
	if (f > 0.0) {
	  RbtDouble phi = 180.0 - Rbt::Angle(cPerp,cAtom1_1,cAtom1_2);
	  RbtDouble Dphi = phi - PHI1prms.R0;
	  f *= f1(fabs(Dphi),PHI1prms);
	  //cout << "LP1:\t" << theta << "," << PHI1prms.R0 << "," << Dphi << ", f=" << f << endl;
	}
      }
      if (f > 0.0) {
	//For guanidinium bPlane1 interacting with C=O lone pair bLP2, we want to use the regular angular dependence
	if (bAngle2 || (bPlane1 && bLP2)) {
	  const RbtCoord& cAtom2_2 = pAtom2_2->GetCoords();
	  RbtDouble DA2 = Rbt::Angle(cAtom1_1,cAtom2_1,cAtom2_2)-A2prms.R0;
	  f *= f1(fabs(DA2),A2prms);
	  //cout << "A2:\t" << A2prms.R0 << "," << DA2 << ", f=" << f << endl;
	}
	else if (bPlane2) {
	  const RbtCoord& cAtom2_2 = pAtom2_2->GetCoords();
	  const RbtCoord& cAtom2_3 = pAtom2_3->GetCoords();
	  RbtPlane pl2 = RbtPlane(cAtom2_1,cAtom2_2,cAtom2_3);
	  RbtDouble A = acos(-fabs(Rbt::Dot(v12.Unit(),pl2.VNorm())))*180.0/M_PI;
	  RbtDouble DA2 = A-A2prms.R0;
	  f *= f1(fabs(DA2),A2prms);
	  //cout << "Pl2:\t" << A2prms.R0 << "," << DA2 << ", f=" << f << endl;
	}
	else if (bLP2) {
	  const RbtCoord& cAtom2_2 = pAtom2_2->GetCoords();
	  const RbtCoord& cAtom2_3 = pAtom2_3->GetCoords();
	  const f1prms& PHI2prms = (eLP2 == RbtInteractionCenter::LONEPAIR) ? m_PHI_lp_prms : m_PHI_plane_prms;
	  RbtPlane pl2 = RbtPlane(cAtom2_1,cAtom2_2,cAtom2_3);
	  RbtDouble dPerp = Rbt::DistanceFromPointToPlane(cAtom1_1, pl2);
	  RbtCoord cPerp = cAtom1_1 - dPerp * pl2.VNorm();
	  //cout << "dPerp = " << dPerp << " check cPerp = " << Rbt::DistanceFromPointToPlane(cPerp, pl2) << endl;
	  RbtDouble theta = asin(dPerp/R)*180.0/M_PI;
	  f *= f1(fabs(theta),m_THETAprms);
	  if (f > 0.0) {
	    RbtDouble phi = 180.0 - Rbt::Angle(cPerp,cAtom2_1,cAtom2_2);
	    RbtDouble Dphi = phi - PHI2prms.R0;
	    f *= f1(fabs(Dphi),PHI2prms);
	    //cout << "LP2:\t" << theta << "," << PHI2prms.R0 << "," << Dphi << ", f=" << f << endl;
	  }
	}
	if (f > 0.0) {
	  //cout << pAtom1_1->GetFullAtomName() << "\t" << pAtom2_1->GetFullAtomName() << "\t" << pAtom2_1->GetUser1Value()
	  //     << "\t" << f << endl;
	  s += pAtom2_1->GetUser1Value()*f;
	  if (bAnnotate) {
	    //DM 6 Feb 2003. Include the charge and receptor density scaling factors in the polar annotation score
	    //so that the score truly reflects the strength of the interaction, as scored by rDock
	    //This also distinguishes attractive (-ve) and repulsive (+ve) annotations
	    RbtAnnotationPtr spAnnotation(new RbtAnnotation(pAtom1_1,pAtom2_1,R,
							    f*pAtom1_1->GetUser1Value()*pAtom2_1->GetUser1Value()));
	    AddAnnotation(spAnnotation);
	  }
	} 
      }
    }
  }
  return s;
}

//As this has a virtual base class we need a separate OwnParameterUpdated
//which can be called by concrete subclass ParameterUpdated methods
//See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
void RbtPolarSF::OwnParameterUpdated(const RbtString& strName) {
  //DM 25 Oct 2000 - heavily used params
  if (strName == _R12FACTOR) {
    m_R12Factor = GetParameter(_R12FACTOR);
  }
  else if (strName == _R12INCR) {
    m_R12Incr = GetParameter(_R12INCR);
  }
  else if (strName == _DR12MIN) {
    m_DR12Min = GetParameter(_DR12MIN);
  }
  else if (strName == _DR12MAX) {
    m_DR12Max = GetParameter(_DR12MAX);
  }
  else if (strName == _A1) {
    m_A1 = GetParameter(_A1);
  }
  else if (strName == _DA1MIN) {
    m_DA1Min = GetParameter(_DA1MIN);
  }
  else if (strName == _DA1MAX) {
    m_DA1Max = GetParameter(_DA1MAX);
  }
  else if (strName == _A2) {
    m_A2 = GetParameter(_A2);
  }
  else if (strName == _DA2MIN) {
    m_DA2Min = GetParameter(_DA2MIN);
  }
  else if (strName == _DA2MAX) {
    m_DA2Max = GetParameter(_DA2MAX);
  }
  else if (strName == _ABS_DR12) {
    m_bAbsDR12 = GetParameter(_ABS_DR12);
  }
  else if (strName == _LP_PHI) {
    m_LP_PHI = GetParameter(_LP_PHI);
    UpdateLPprms();
  }
  else if (strName == _LP_DPHIMIN) {
    m_LP_DPHIMin = GetParameter(_LP_DPHIMIN);
    UpdateLPprms();
  }
  else if (strName == _LP_DPHIMAX) {
    m_LP_DPHIMax = GetParameter(_LP_DPHIMAX);
    UpdateLPprms();
  }
  else if (strName == _LP_DTHETAMIN) {
    m_LP_DTHETAMin = GetParameter(_LP_DTHETAMIN);
    UpdateLPprms();
  }
  else if (strName == _LP_DTHETAMAX) {
    m_LP_DTHETAMax = GetParameter(_LP_DTHETAMAX);
    UpdateLPprms();
  }
}


void RbtPolarSF::UpdateLPprms() {
  m_PHI_lp_prms = f1prms(m_LP_PHI,m_LP_DPHIMin,m_LP_DPHIMax);
  m_PHI_plane_prms = f1prms(0.0,m_LP_PHI+m_LP_DPHIMin,m_LP_PHI+m_LP_DPHIMax);
  m_THETAprms = f1prms(0.0,m_LP_DTHETAMin,m_LP_DTHETAMax);
}

