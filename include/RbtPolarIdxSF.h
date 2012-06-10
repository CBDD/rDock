/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPolarIdxSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Indexed-grid-based intermolecular scoring function for all attractive polar
//interactions (HBD,HBA,metal,guanidinium carbon)

#ifndef _RBTPOLARIDXSF_H_
#define _RBTPOLARIDXSF_H_

#include "RbtBaseIdxSF.h"
#include "RbtBaseInterSF.h"
#include "RbtPolarSF.h"

class RbtPolarIdxSF : public RbtBaseInterSF, public RbtBaseIdxSF, public RbtPolarSF
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _INCR;
  static RbtString _ATTR;
  //DM 12 Jun 2002 - score thresholds used for counting attractive and repulsive interactions
  static RbtString _THRESHOLD_POS;
  static RbtString _THRESHOLD_NEG;
  
  RbtPolarIdxSF(const RbtString& strName = "POLAR");
  virtual ~RbtPolarIdxSF();

  //Override RbtBaseSF::ScoreMap to provide additional raw descriptors
  virtual void ScoreMap(RbtStringVariantMap& scoreMap) const;  

 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupSolvent();
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  
  //Clear the receptor and ligand grids and lists respectively
  //As we are not using smart pointers, there is some memory management to do
  void ClearReceptor();
  void ClearLigand();
  void ClearSolvent();
  //Helper function for above
  void DeleteList(RbtInteractionCenterList& icList);
  
  //DM 25 Oct 2000 - track changes to parameter values in local data members
  //ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const RbtString& strName);
  
 private:
  RbtDouble ReceptorScore() const;
  RbtDouble SolventScore() const;
  RbtDouble InterScore() const;
  RbtDouble ReceptorSolventScore() const;
  RbtDouble LigandSolventScore() const;
  
  RbtDouble InterScore(const RbtInteractionCenterList& posList,
  						const RbtInteractionCenterList& negList,
  						RbtBool bCount) const;
  RbtInteractionGridPtr m_spPosGrid;
  RbtInteractionGridPtr m_spNegGrid;
  RbtInteractionCenterList m_recepPosList;
  RbtInteractionCenterList m_recepNegList;
  RbtInteractionCenterList m_flexRecPosList;
  RbtInteractionCenterList m_flexRecNegList;
  RbtInteractionListMap m_flexRecIntns;
  RbtInteractionListMap m_flexRecPrtIntns;//partitioned interactions

  RbtInteractionCenterList m_ligPosList;
  RbtInteractionCenterList m_ligNegList;
  
  RbtInteractionCenterList m_solventPosList;
  RbtInteractionCenterList m_solventNegList;
  RbtInteractionListMap m_solventIntns;

  RbtBool m_bAttr;
  RbtBool m_bFlexRec;
  RbtBool m_bSolvent;
  //DM 12 Jun 2002 - keep track of number of ligand atoms involved in non-zero polar interactions
  mutable RbtInt m_nPos;//#positive centers with non-zero scores
  mutable RbtInt m_nNeg;//#negative centers with non-zero scores
  RbtDouble m_posThreshold;
  RbtDouble m_negThreshold;
};

#endif //_RBTPOLARIDXSF_H_
