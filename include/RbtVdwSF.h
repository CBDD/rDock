/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtVdwSF.h#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Base implementation class for all vdW scoring functions

#ifndef _RBTVDWSF_H_
#define _RBTVDWSF_H_

#include "RbtBaseSF.h"
#include "RbtTriposAtomType.h"
#include "RbtAtom.h"
#include "RbtParameterFileSource.h"
#include "RbtAnnotationHandler.h"

class RbtVdwSF : public virtual RbtBaseSF, public virtual RbtAnnotationHandler
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _USE_4_8;//TRUE = 4-8; FALSE = 6-12
  static RbtString _USE_TRIPOS;//TRUE = Tripos 5.2 well depths; FALSE = GOLD well depths
  static RbtString _RMAX;//Maximum range as a multiple of rmin
  static RbtString _ECUT;//Energy cutoff for transition to short-range quadratic, as a multiple of well depth
  static RbtString _E0;//Energy at zero distance, as a multiple of ECUT
  
  virtual ~RbtVdwSF();
  
 protected:
  RbtVdwSF();
  
  //As this has a virtual base class we need a separate OwnParameterUpdated
  //which can be called by concrete subclass ParameterUpdated methods
  //See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
  void OwnParameterUpdated(const RbtString& strName);

  //Used by subclasses to calculate vdW potential between pAtom and all atoms in atomList
  RbtDouble VdwScore(const RbtAtom* pAtom, const RbtAtomRList& atomList) const;
  //As above, but with additional checks for enabled state of each atom
  RbtDouble VdwScoreEnabledOnly(const RbtAtom* pAtom, const RbtAtomRList& atomList) const;
  //Looks up the maximum range (rmax_sq) for any interaction
  //i.e. from across a row of m_vdwTable
  RbtDouble MaxVdwRange(const RbtAtom* pAtom) const;
  RbtDouble MaxVdwRange(RbtTriposAtomType::eType t) const;

  //Index the intramolecular interactions
  void BuildIntraMap(const RbtAtomRList& atomList, RbtAtomRListList& intns) const;
  void BuildIntraMap(const RbtAtomRList& atomList1, const RbtAtomRList& atomList2, RbtAtomRListList& intns) const;
  void Partition(const RbtAtomRList& atomList, const RbtAtomRListList& intns, RbtAtomRListList& prtIntns, RbtDouble dist=0.0) const;

 private:
  //vdW scoring function params
  struct vdwprms {
    RbtDouble A,B;//6-12 or 4-8 params
    RbtDouble rmin;//Sum of vdw radii
    RbtDouble kij;//Well depth
    RbtDouble rmax_sq;//Max distance**2 over which the vdw score should be calculated
    RbtDouble rcutoff_sq;//Distance**2 at which the short-range quadratic potential kicks in
    RbtDouble ecutoff;//Energy at the cutoff point
    RbtDouble slope;//Slope of the short-range quadratic potential
    RbtDouble e0;//Energy at zero distance
  };

  typedef vector<RbtVdwSF::vdwprms> RbtVdwRow;
  typedef RbtVdwRow::iterator RbtVdwRowIter;
  typedef RbtVdwRow::const_iterator RbtVdwRowConstIter;  
  typedef vector<RbtVdwRow> RbtVdwTable;
  typedef RbtVdwTable::iterator RbtVdwTableIter;
  typedef RbtVdwTable::const_iterator RbtVdwTableConstIter;

  //Generic scoring function primitive for 6-12
  inline RbtDouble f6_12(RbtDouble R_sq, const vdwprms& prms) const {
    //Zero well depth or long range: return zero
    if ( (prms.kij == 0.0) || (R_sq > prms.rmax_sq) ) {
      return 0.0;
    }
    //Short range: use quadratic potential
    else if (R_sq < prms.rcutoff_sq) {
      return prms.e0 - (prms.slope * R_sq);
    }
    //Everywhere else is 6-12
    else {
      RbtDouble rr6 = 1.0 / (R_sq * R_sq * R_sq);
      return rr6 * (rr6 * prms.A - prms.B);
    }
  };
  
  //Generic scoring function primitive for 4-8
  inline RbtDouble f4_8(RbtDouble R_sq, const vdwprms& prms) const {
    //Zero well depth or long range: return zero
    if ( (prms.kij == 0.0) || (R_sq > prms.rmax_sq) ) {
      return 0.0;
    }
    //Short range: use quadratic potential
    else if (R_sq < prms.rcutoff_sq) {
      return prms.e0 - (prms.slope * R_sq);
    }
    //Everywhere else is 4-8
    else {
      RbtDouble rr4 = 1.0 / (R_sq * R_sq);
      return rr4 * (rr4 * prms.A - prms.B);
    }
  };

  void Setup();//Initialise m_vdwTable with appropriate params for each atom type pair
  void SetupCloseRange();//Regenerate the short-range params only (called more frequently)

  //Private predicate
  //Is the distance between atoms less than a given value ?
  //Function checks d**2 to save performing a sqrt
  class isD_lt : public std::unary_function<RbtAtom*,RbtBool> {
    RbtDouble d_sq;
    RbtAtom* a;
  public:
    explicit isD_lt(RbtAtom* aa, RbtDouble dd) : d_sq(dd*dd),a(aa) {}
    RbtBool operator() (RbtAtom* aa) const {return Rbt::Length2(aa->GetCoords(),a->GetCoords()) < d_sq;}
  };

  //Private data members
  RbtBool m_use_4_8;
  RbtBool m_use_tripos;
  RbtDouble m_rmax;
  RbtDouble m_ecut;
  RbtDouble m_e0;
  RbtParameterFileSourcePtr m_spVdwSource;//File source for vdw params
  RbtVdwTable m_vdwTable;//Lookup table for all vdW params (indexed by Tripos atom type)
  RbtDoubleList m_maxRange;//Vector of max ranges for each Tripos atom type
};

#endif //_RBTVDWSF_H_
