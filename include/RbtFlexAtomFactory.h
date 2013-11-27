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

//RbtFlexDataVisitor class to generate lists of fixed, tethered and free atoms
//from a model or list of models
#ifndef _RBTFLEXATOMFACTORY_H_
#define _RBTFLEXATOMFACTORY_H_

#include "RbtFlexDataVisitor.h"
#include "RbtModel.h"

//The main use of RbtFlexAtomFactory is to aid with scoring function performance
//enhancements by partitioning model atoms into subsets that are:
// - fixed (never move)
// - tethered (move, but within a predictable radius)
// - free (move freely)
//
//For receptor models, the majority of atoms are FIXED
//The flexible OH/NH3+ terminal atoms are defined as TETHERED
//
//Two categorisations are defined for ligands and solvent 
//depending on the internal flexibility of the model
//
//A) If model is flexible (RbtModel::isFlexible()), all atoms are defined as FREE
//as it is too difficult to take the impact of internal dihedral motion into account
//
//B) If model is not flexible (e.g. explicit solvent), the mapping is:
//transMode == FIXED && rotMode == FIXED => atoms are FIXED
//transMode == FREE => atoms are FREE
//else atoms are TETHERED
//
//On exit:
//User2Value of each tethered atom is set to the maximum displacement
//possible for translational and rotational moves. The value is a
//conservative overestimate as the allowed volume mapped on by each tethered
//atom can be complex. Here we define the sphere radius that includes all
//allowed coordinates.
//The User2Value for fixed and free atoms is set to zero.
class RbtFlexAtomFactory : public RbtFlexDataVisitor {
 public:
  RbtFlexAtomFactory() {};
  //Constructor accepting a single model
  RbtFlexAtomFactory(RbtModel*);
  //Constructor accepting a list of models
  RbtFlexAtomFactory(RbtModelList);
  //Helper method to visit a single model
  void Visit(RbtModel*);
  //Clears the partitioned atom lists
  void Clear();
  //Implement the base class abstract methods
  virtual void VisitReceptorFlexData(RbtReceptorFlexData*);
  virtual void VisitLigandFlexData(RbtLigandFlexData*);
  virtual void VisitSolventFlexData(RbtSolventFlexData*);
  //Get methods for the partitioned atom lists
  RbtAtomRList GetFixedAtomList() const {return m_fixedAtomList;}
  RbtAtomRList GetTetheredAtomList() const {return m_tetheredAtomList;}
  RbtAtomRList GetFreeAtomList() const {return m_freeAtomList;}
 private:
  RbtAtomRList m_fixedAtomList;
  RbtAtomRList m_tetheredAtomList;
  RbtAtomRList m_freeAtomList;
};
#endif //_RBTFLEXATOMFACTORY_H_
