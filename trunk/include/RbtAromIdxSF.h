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

//Indexed-grid-based intermolecular aromatic scoring function

#ifndef _RBTAROMIDXSF_H_
#define _RBTAROMIDXSF_H_

#include "RbtBaseIdxSF.h"
#include "RbtBaseInterSF.h"
#include "RbtPlane.h"
#include "RbtAnnotationHandler.h"

class RbtAromIdxSF : public RbtBaseInterSF, public RbtBaseIdxSF, public RbtAnnotationHandler
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _INCR;
  static RbtString _R12;
  static RbtString _DR12MIN;
  static RbtString _DR12MAX;
  static RbtString _DAMIN;
  static RbtString _DAMAX;
  //DM 12 Jun 2002 - score threshold used for counting aromatic interactions
  static RbtString _THRESHOLD;
		
  RbtAromIdxSF(const RbtString& strName = "AROM");
  virtual ~RbtAromIdxSF();

  //Override RbtBaseSF::ScoreMap to provide additional raw descriptors
  //virtual void ScoreMap(RbtStringVariantMap& scoreMap) const;  
  
 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;

  //Clear the receptor and ligand grids and lists respectively
  //As we are not using smart pointers, there is some memory management to do
  void ClearReceptor();
  void ClearLigand();

  //DM 25 Oct 2000 - track changes to parameter values in local data members
  //ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const RbtString& strName);
  
 private:
  //////////////////////////////////////////////////////////
  //DM 6 Feb 2003
  //To be consistent with the Polar score, we should really define a base RbtAromSF class
  //which will provide the scoring function primitives, then have an RbtAromIdxSF subclass
  //for intermolecular aromatic interactions, and a RbtAromIntraSF subclass for
  //ligand intramolecular interactions. Currently there is no RbtAromIntraSF class so
  //everything is contained within RbtAromIdxSF

  //Generic scoring function params
  struct f1prms {
    RbtDouble R0,DRMin,DRMax,slope;
    f1prms::f1prms(RbtDouble R, RbtDouble DMin, RbtDouble DMax)
      : R0(R),DRMin(DMin),DRMax(DMax),slope(1.0/(DMax-DMin)) {};
  };

  inline f1prms GetRprms() const {return f1prms(m_R12,m_DR12Min,m_DR12Max);}
  inline f1prms GetAprms() const {return f1prms(0.0,m_DAMin,m_DAMax);}

  //Generic scoring function primitive (deviation from ideal geometry)
  inline RbtDouble f1(RbtDouble DR, const f1prms& prms) const {
    return (DR > prms.DRMax) ? 0.0 : (DR > prms.DRMin) ? 1.0-prms.slope*(DR-prms.DRMin) : 1.0;
  };

  //The actual aromatic score, between a given interaction center and a list of near neighbour centers
  RbtDouble AromScore(const RbtInteractionCenter* pIC1, const RbtInteractionCenterList& IC2List,
		      const f1prms& Rprms, const f1prms& Aprms) const;
  RbtDouble PiScore(const RbtInteractionCenter* pIC1, const RbtInteractionCenterList& IC2List) const;
  //End of section that should ultimately be moved to RbtAromSF base class
  //////////////////////////////////////////////////////////

  RbtInteractionGridPtr m_spAromGrid;
  RbtInteractionGridPtr m_spGuanGrid;
  RbtInteractionCenterList m_recepAromList;
  RbtInteractionCenterList m_recepGuanList;
  RbtInteractionCenterList m_ligAromList;
  RbtInteractionCenterList m_ligGuanList;
  
  //DM 25 Oct 2000 - local copies of params
  RbtDouble m_R12;
  RbtDouble m_DR12Min;
  RbtDouble m_DR12Max;
  RbtDouble m_DAMin;
  RbtDouble m_DAMax;

  //DM 12 Jun 2002 - keep track of number of ligand rings and guan carbons involved in non-zero arom interactions
  mutable RbtInt m_nArom;
  mutable RbtInt m_nGuan;
  mutable RbtDouble m_ss;//Sigma-sigma score
  mutable RbtDouble m_sp;//Sigma-pi score
  mutable RbtDouble m_pp;//pi-pi score 
  RbtDouble m_threshold;
};

#endif //_RBTAROMIDXSF_H_
