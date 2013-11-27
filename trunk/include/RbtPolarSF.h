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

//Base implementation class for all attractive and repulsive polar scoring fns
//Provides methods for: creating interaction centers, calculating interaction
//scores

#ifndef _RBTPOLARSF_H_
#define _RBTPOLARSF_H_

#include "RbtBaseSF.h"
#include "RbtInteractionGrid.h"
#include "RbtAnnotationHandler.h"

class RbtPolarSF : public virtual RbtBaseSF, public virtual RbtAnnotationHandler
{
 public:
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _INCR;
  static RbtString _R12FACTOR;
  static RbtString _R12INCR;
  static RbtString _DR12MIN;
  static RbtString _DR12MAX;
  static RbtString _A1;
  static RbtString _DA1MIN;
  static RbtString _DA1MAX;
  static RbtString _A2;
  static RbtString _DA2MIN;
  static RbtString _DA2MAX;  
  static RbtString _INCMETAL;
  static RbtString _INCHBD;
  static RbtString _INCHBA;
  static RbtString _INCGUAN;
  static RbtString _GUAN_PLANE;
  static RbtString _ABS_DR12;
  static RbtString _LP_OSP2;
  static RbtString _LP_PHI;
  static RbtString _LP_DPHIMIN;
  static RbtString _LP_DPHIMAX;
  static RbtString _LP_DTHETAMIN;
  static RbtString _LP_DTHETAMAX;
  
  virtual ~RbtPolarSF();
  
 protected:
  RbtPolarSF();

  RbtInteractionCenterList CreateAcceptorInteractionCenters(const RbtAtomList& atomList) const;
  RbtInteractionCenterList CreateDonorInteractionCenters(const RbtAtomList& atomList) const;

  //Index the intramolecular interactions between two lists
  void BuildIntraMap(const RbtInteractionCenterList& ICList1,
		     const RbtInteractionCenterList& ICList2,
		     RbtInteractionListMap& intns) const;
  //Index the intramolecular interactions within a single list
  void BuildIntraMap(const RbtInteractionCenterList& ICList,
		     RbtInteractionListMap& intns) const;

  RbtDouble IntraScore(const RbtInteractionCenterList& posList,
		       const RbtInteractionCenterList& negList,
		       const RbtInteractionListMap& prtIntns, RbtBool attr) const;
  void Partition(const RbtInteractionCenterList& posList,
		 const RbtInteractionCenterList& negList,
		 const RbtInteractionListMap& intns,
		 RbtInteractionListMap& prtIntns, RbtDouble dist=0.0) const;
  
  //Generic scoring function params
  struct f1prms {
    RbtDouble R0,DRMin,DRMax,slope;
    f1prms::f1prms(RbtDouble R, RbtDouble DMin, RbtDouble DMax)
      : R0(R),DRMin(DMin),DRMax(DMax),slope(1.0/(DMax-DMin)) {};
  };

  inline f1prms GetRprms() const {return f1prms(0.0,m_DR12Min,m_DR12Max);}
  inline f1prms GetA1prms() const {return f1prms(m_A1,m_DA1Min,m_DA1Max);}
  inline f1prms GetA2prms() const {return f1prms(m_A2,m_DA2Min,m_DA2Max);}

  RbtDouble PolarScore(const RbtInteractionCenter* intn, const RbtInteractionCenterList& intnList,
		       const f1prms& Rprms, const f1prms& A1prms, const f1prms& A2prms) const;

  //As this has a virtual base class we need a separate OwnParameterUpdated
  //which can be called by concrete subclass ParameterUpdated methods
  //See Stroustrup C++ 3rd edition, p395, on programming virtual base classes
  void OwnParameterUpdated(const RbtString& strName);

 private:
  //Generic scoring function primitive
  inline RbtDouble f1(RbtDouble DR, const f1prms& prms) const {
    return (DR > prms.DRMax) ? 0.0 : (DR > prms.DRMin) ? 1.0-prms.slope*(DR-prms.DRMin) : 1.0;
  };

  void UpdateLPprms();

  //DM 25 Oct 2000 - heavily used params
  RbtDouble m_R12Factor;
  RbtDouble m_R12Incr;
  RbtDouble m_DR12Min;
  RbtDouble m_DR12Max;
  RbtDouble m_A1;
  RbtDouble m_DA1Min;
  RbtDouble m_DA1Max;
  RbtDouble m_A2;
  RbtDouble m_DA2Min;
  RbtDouble m_DA2Max;
  RbtBool m_bAbsDR12;
  RbtDouble m_LP_PHI;
  RbtDouble m_LP_DPHIMin;
  RbtDouble m_LP_DPHIMax;
  RbtDouble m_LP_DTHETAMin;
  RbtDouble m_LP_DTHETAMax;
  f1prms m_PHI_lp_prms;
  f1prms m_PHI_plane_prms;
  f1prms m_THETAprms;
};

#endif //_RBTPOLARSF_H_
