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

//Desolvation scoring function based on Weighted Solvent Accessible Surface Area

#ifndef _RBTSAIDXSF_H_
#define _RBTSAIDXSF_H_

#include "RbtBaseInterSF.h"
#include "RbtBaseIdxSF.h"
#include "RbtSATypes.h"
#include "RbtNonBondedHHSGrid.h"
#include "RbtParameterFileSource.h"
#include "RbtAnnotationHandler.h"

class RbtSAIdxSF : public RbtBaseInterSF, public RbtBaseIdxSF, public RbtAnnotationHandler
{
 public:
  RbtSAIdxSF(const RbtString& strName="SAIdxSF");
  virtual ~RbtSAIdxSF();
  // write score components
  virtual void ScoreMap(RbtStringVariantMap& scoreMap) const;
  
  static RbtString _CT;
  static RbtString _INCR;

  //Request Handling method
  //Handles the Partition request
  virtual void HandleRequest(RbtRequestPtr spRequest);

 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupSolvent();
  virtual void SetupScore();
  virtual RbtDouble RawScore(void) const;
  
  // clear indexed grids
  void ClearReceptor(void);
  void ClearLigand(void);
  void ClearSolvent(void);
  
  RbtDouble GetASP(RbtHHSType::eType, RbtDouble) const;
  RbtDouble GetP_i(RbtHHSType::eType) const;
  RbtDouble GetR_i(RbtHHSType::eType) const;
  void PrintWeightMatrix() const;

 private:
  struct solvprms {
    solvprms(RbtDouble r0, RbtDouble p0, RbtDouble asp0, RbtBool b) :
      r(r0),p(p0),asp(asp0),chg_scaling(b) {};
    RbtDouble r;
    RbtDouble p;
    RbtDouble asp;
    RbtBool chg_scaling;
  };

  typedef vector<RbtSAIdxSF::solvprms> RbtSolvTable;
  void Setup();
  HHS_SolvationRList CreateInteractionCenters(const RbtAtomList& atomList) const;
  void BuildIntraMap(HHS_SolvationRList& ICList) const;
  void BuildIntraMap(HHS_SolvationRList& ICList1, HHS_SolvationRList& ICList2) const;

  //Sum the surface energies (ASP*area) for the list of solvation interaction centers
  RbtDouble TotalEnergy(const HHS_SolvationRList& intnCenters) const;
  void Partition(HHS_SolvationRList& intnCenters, RbtDouble dist=0.0);

  HHS_SolvationRList theLSPList;//All ligand solvation interaction centers
  HHS_SolvationRList theRSPList;//All rigid receptor solvation interaction centers
  HHS_SolvationRList theFlexList;//All flexible receptor solvation interaction centers
  HHS_SolvationRList theCavList;//Subset of rigid receptor list within range of docking site
  //Subset of rigid receptor list NOT within range of docking site
  //but within range of the flexible atoms
  //Must remember to restore the invariant areas for these atoms
  //before updating the overlap of the flexible atoms
  HHS_SolvationRList thePeriphList;
  HHS_SolvationRList theSolventList;//DM 21 Dec 2005 - explicit solvent interaction centers
  RbtNonBondedHHSGridPtr theIdxGrid;
  RbtSolvTable m_solvTable;
  RbtParameterFileSourcePtr m_spSolvSource;//File source for solvation params
  RbtDouble m_maxR;//Maximum radius of any atom type, used to adjust Range() dynamically
  RbtBool m_bFlexRec;//Is receptor flexible?
  mutable RbtDouble m_lig_0;//Solvation energy of the free ligand (initial conformation)
  mutable RbtDouble m_lig_free;//Solvation energy of the free ligand (current conformation)
  mutable RbtDouble m_lig_bound;//Solvation energy of the bound ligand (current conformation)
  mutable RbtDouble m_site_0;//Solvation energy of the free docking site (initial conformation)
  mutable RbtDouble m_site_free;//Solvation energy of the free docking site (current conformation)
  mutable RbtDouble m_site_bound;//Solvation energy of the bound docking site (current conformation)
  mutable RbtDouble m_solvent_0;//Solvation energy of the free explicit solvent (initial conformation)
  mutable RbtDouble m_solvent_free;//Solvation energy of the free explicit solvent (current conformation)
  mutable RbtDouble m_solvent_bound;//Solvation energy of the bound explicit solvent (current conformation)
};

#endif // _RBTSAIDXSF_H_ 
