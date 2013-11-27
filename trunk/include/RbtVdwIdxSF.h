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

//Indexed-grid-based intermolecular scoring function for vdw interactions
#ifndef _RBTVDWIDXSF_H_
#define _RBTVDWIDXSF_H_

#include "RbtBaseIdxSF.h"
#include "RbtBaseInterSF.h"
#include "RbtVdwSF.h"

class RbtVdwIdxSF : public RbtBaseInterSF, public RbtBaseIdxSF, public RbtVdwSF
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  //DM 12 Jun 2002 - score thresholds used for counting attractive and repulsive interactions
  //(thresholds represent partial vdW scores summed per ligand atom,
  //used for outputting ligand atom counts)
  static RbtString _THRESHOLD_ATTR;
  static RbtString _THRESHOLD_REP;
  //DM 03 Feb 2003 - score threshold for outputting individual (atom-atom) lipophilic vdW interactions
  //i.e between two non-polar carbons/hydrogens
  static RbtString _ANNOTATION_LIPO;
  //DM 10 Apr 2003 - vdw annotation can be quite slow so provide boolean option to disable
  static RbtString _ANNOTATE;
  //DM 14 Jun 2006 - option to disable solvent performance enhancements, mainly for testing
  static RbtString _FAST_SOLVENT;
  
  RbtVdwIdxSF(const RbtString& strName = "VDW");
  virtual ~RbtVdwIdxSF();
  
  //Override RbtBaseSF::ScoreMap to provide additional raw descriptors
  virtual void ScoreMap(RbtStringVariantMap& scoreMap) const;

 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupSolvent();
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  RbtDouble InterScore() const;
  RbtDouble ReceptorScore() const;
  RbtDouble SolventScore() const;
  RbtDouble ReceptorSolventScore() const;
  RbtDouble LigandSolventScore() const;

  //DM 25 Oct 2000 - track changes to parameter values in local data members
  //ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const RbtString& strName);
  
 private:
  void RenderAnnotationsByResidue(RbtStringList& retVal) const;

  RbtNonBondedGridPtr m_spGrid;//Indexing grid for receptor
  RbtNonBondedGridPtr m_spSolventGrid;//Indexing grid for fixed/tethered solvent
  RbtAtomList m_recAtomList;
  RbtAtomRList m_recRigidAtomList;
  RbtAtomRList m_recFlexAtomList;
  RbtAtomRList m_ligAtomList;
  RbtAtomRList m_solventAtomList;//All solvent atoms
  RbtAtomRList m_solventFixTethAtomList;//Fixed or tethered solvent atoms
  RbtAtomRList m_solventFreeAtomList;//Free solvent atoms
  RbtAtomRListList m_recFlexIntns;//Intra-target flexible interactions
  RbtAtomRListList m_recFlexPrtIntns;//Intra-target flexible partitioned interactions
  RbtAtomRListList m_solventFixTethIntns;//Intra-solvent intns between fixed/tethered atoms
  RbtAtomRListList m_solventFixTethPrtIntns;//Partitioned intns between fixed/tethered solvent
  RbtAtomRListList m_solventFreeIntns;//Intra-solvent intns between free solvent atoms
  //DM 12 Jun 2002 - keep track of number of ligand atoms involved in non-zero vdW interactions
  mutable RbtInt m_nAttr;//#atoms with net attractive (-ve) vdw scores
  mutable RbtInt m_nRep;//#atoms with net repulsive (+ve) vdw scores
  RbtDouble m_attrThreshold;
  RbtDouble m_repThreshold;
  RbtDouble m_lipoAnnot;
  RbtBool m_bAnnotate;
  RbtBool m_bFlexRec;
  RbtBool m_bFastSolvent;
};

#endif //_RBTVDWIDXSF_H_
