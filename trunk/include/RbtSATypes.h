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

#ifndef _RBTSATYPES_H_
#define _RBTSATYPES_H_

#include "RbtConfig.h"

class RbtAtom;

//Atom typing class for the solvation atom types
//HHS = Hasel, Hendrickson and Still, Tetrahedron Computer Methodology, 1 (1988) 103-116
//This paper describes the fast approximation to solvent-accessible surface area used here
//Solvation atom types are inspired by Wang, Wang, Huo, Lee and Kollman, J. Phys. Chem. B 105 (2001) 5055-5067
class RbtHHSType {
 public:
  enum eType {
    UNDEFINED = 0,
    C_sp3, CH_sp3, CH2_sp3, CH3_sp3,
    C_sp2, CH_sp2, CH2_sp2, C_sp2p,
    C_ar, CH_ar,
    C_sp,
    C_sp3_P, CH_sp3_P, CH2_sp3_P, CH3_sp3_P,
    C_sp2_P, CH_sp2_P, CH2_sp2_P,
    C_ar_P, CH_ar_P,
    H, HO, HN, HNp, HS,
    O_sp3, OH_sp3, OW, O_tri, OH_tri, O_sp2, ON, Om,
    N_sp3, NH_sp3, NH2_sp3, N_sp3p, N_tri, NH_tri, NH2_tri, N_sp2, N_sp2p, N_ar, N_sp,
    S_sp3, S_sp2,
    P,
    F, Cl, Br, I,
    Metal,
    MAXTYPES //KEEP AS LAST TYPE: used to size the atom name string list
  };
  
  RbtHHSType();
  virtual ~RbtHHSType();
  
  eType	operator() (RbtAtom* anAtom) const;
  RbtString Type2Str(eType) const;
  eType	Str2Type(const RbtString&) const;

 private:  
  void SetupTypeNames();
  RbtStringList m_typeNames;
};

typedef SmartPtr<RbtHHSType*> RbtHHSTypePtr;	//Smart pointer
typedef vector<RbtHHSType::eType> RbtHHSTypeList;
typedef RbtHHSTypeList::iterator RbtHHSTypeListIter;
typedef RbtHHSTypeList::const_iterator RbtHHSTypeListConstIter;

//Interaction center class for solvation calculations
//A "decorated" RbtAtom in effect
//Designed as a "create-once, read-many" object
//Attributes (radius, ASP etc) cannot be changed once constructed
//The Overlap method accumulates the surface areas for each atom
//based on pair-wise distances to surrounding atoms
class HHS_Solvation 
{
 public:
  //Constants from the HHS paper
  static const RbtDouble r_s;//solvent probe radius
  static const RbtDouble d_s;//solvent probe diameter
  static const RbtDouble Pij_12; // correction factor for 1-2 connections
  static const RbtDouble Pij_13; // ditto for 1-3
  static const RbtDouble Pij_14; // and for 1-4+

  //Default p_i, r_i, and sigma values
  HHS_Solvation(RbtHHSType::eType t, RbtAtom* a, RbtDouble p=1.0, RbtDouble r=1.0, RbtDouble s=0.0);

  //Get properties
  RbtDouble GetA_i(void) const {return A_i;}
  RbtDouble GetS_i(void) const {return S_i;}
  RbtDouble GetP_i(void) const {return p_i;}
  RbtDouble GetR_i(void) const {return r_i;}
  RbtDouble GetSigma(void) const {return sigma;}
  RbtHHSType::eType GetHHSType() const {return hhsType;};
  RbtAtom* GetAtom() const {return atom;};
  inline RbtDouble GetArea() const {return S_i * A_i;};//Exposed area
  inline RbtDouble GetEnergy() const {return E_i * A_i;};//Surface energy

  //Actions
  //Initialise all calculated parameters to an isolated atom (no overlap)
  void Init();
  //Saves current value of A_i to A_inv
  //Use to memorise the partially calculated area due to invariant interaction distances
  inline void Save() {A_inv = A_i;};
  //Restores A_inv to A_i
  //Use prior to continuing the calculation due to variable interaction distances
  inline void Restore() {A_i = A_inv;};
  //Calculate overlap between this center and another (h)
  //Updates the exposed fractions (A_i) for both centers
  //p_ij is the correction factor for 1-2, 1-3, and 1-4+ connected atoms
  void Overlap(HHS_Solvation* h, RbtDouble p_ij);
  //Get a list of all the current (partitioned) variable-distance interactions to this center
  const vector<HHS_Solvation*>& GetVariable() const {return m_prt;}
  RbtInt GetNumVariable() const {return m_prt.size();}
  //Add a variable-distance interaction
  void AddVariable(HHS_Solvation* anAtom);
  //Helper function to calculate the overlap for all the variable interactions to this center
  //All are 1-4+ by definition
  void OverlapVariable();
  //As above, but ignores disabled interaction centers
  void OverlapVariableEnabledOnly();
  //Partition the list of variable interactions
  //Copy those that are closer than distance d to the partitioned list
  //d=0 => clear the partition
  void Partition(RbtDouble d=0.0);
  
 private:
  RbtDouble p_i;//atom-type dependent correction
  RbtDouble r_i;//quasi vdW radius
  RbtDouble S_i;//Surface area of isolated atom (constant)
  RbtDouble A_i;//Fraction of atom surface currently exposed
  RbtDouble A_inv;//invariant part of atom surface (memory, used by Save/Restore)
  RbtDouble sigma;//atomic solvation parameter
  RbtDouble E_i;//energy of isolated atom (S_i * sigma)
  RbtDouble PI_r_i_plus_r_s;//PI*(r_i+r_s)
  RbtDouble p_i_over_S_i;//p_i / S_i

  RbtHHSType::eType hhsType;//solvation atom type
  RbtAtom* atom;//the RbtAtom itself
  vector<HHS_Solvation*> m_var;//Vector of all variable distances
  vector<HHS_Solvation*> m_prt;//Vector of current partioned variable distances
};

typedef vector<HHS_Solvation*> HHS_SolvationRList;
typedef HHS_SolvationRList::iterator HHS_SolvationRListIter;
typedef HHS_SolvationRList::const_iterator HHS_SolvationRListConstIter;

typedef vector<HHS_SolvationRList> HHS_SolvationListMap;// vector of regular pointers
typedef HHS_SolvationListMap::iterator HHS_SolvationListMapIter;
typedef HHS_SolvationListMap::const_iterator HHS_SolvationListMapConstIter;

typedef SmartPtr<HHS_Solvation> HHS_SolvationPtr;//Smart pointer
typedef vector<HHS_SolvationPtr> HHS_SolvationList;//Vector of smart pointers
typedef HHS_SolvationList::iterator HHS_SolvationListIter;
typedef HHS_SolvationList::const_iterator HHS_SolvationListConstIter;

//Useful STL predicates and function objects
namespace Rbt {
  class InitHHS {
  public:
    explicit InitHHS() {}
    void operator() (HHS_Solvation* pHHS) {pHHS->Init();}
  };
  class SaveHHS {
  public:
    explicit SaveHHS() {}
    void operator() (HHS_Solvation* pHHS) {pHHS->Save();}
  };
  class RestoreHHS {
  public:
    explicit RestoreHHS() {}
    void operator() (HHS_Solvation* pHHS) {pHHS->Restore();}
  };
  class PartitionHHS {
    RbtDouble d;
  public:
    explicit PartitionHHS(RbtDouble dd) : d(dd) {}
    void operator() (HHS_Solvation* pHHS) {pHHS->Partition(d);}
  };
  class OverlapVariableHHS {
  public:
    explicit OverlapVariableHHS() {}
    void operator() (HHS_Solvation* pHHS) {pHHS->OverlapVariable();}
  };
  class OverlapHHS {
    HHS_Solvation* pHHSi;
    RbtDouble p_ij;
  public:
    explicit OverlapHHS(HHS_Solvation* pHHS, RbtDouble p) : pHHSi(pHHS),p_ij(p) {}
    void operator() (HHS_Solvation* pHHSj) {pHHSi->Overlap(pHHSj,p_ij);}
  };

  //Is HHS type equal to t ?
  class isHHSType_eq : public std::unary_function<HHS_Solvation*,RbtBool> {
    RbtHHSType::eType t;
  public:
    explicit isHHSType_eq(RbtHHSType::eType tt) : t(tt) {}
    RbtBool operator() (const HHS_Solvation* pHHS) const {return pHHS->GetHHSType() == t;}
  };

  //Is HHS_Solvation selected ?
  class isHHSSelected : public std::unary_function<HHS_Solvation*,RbtBool> {
  public:
    explicit isHHSSelected() {}
    RbtBool operator() (const HHS_Solvation* pHHS) const;
  };

  inline RbtDouble AccumHHSEnergy(RbtDouble val, const HHS_Solvation* pHHS) {
    return val+pHHS->GetEnergy();
  }
}
#endif //_RBTSATYPES_H_


