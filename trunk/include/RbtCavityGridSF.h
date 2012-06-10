/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCavityGridSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Penalty function for keeping ligand and solvent within the docking volume
#ifndef _RBTCAVITYGRIDSF_H_
#define _RBTCAVITYGRIDSF_H_

#include "RbtBaseInterSF.h"
#include "RbtRealGrid.h"
#include "RbtFlexDataVisitor.h"

class RbtCavityGridSF : public RbtBaseInterSF
{
 public:
    //RbtFlexDataVisitor class to generate list of movable non-hydrogen atoms
    //i.e. those atoms that have the potential to move outside of the cavity
    class HeavyAtomFactory : public RbtFlexDataVisitor {
    public:
        HeavyAtomFactory() {};
        HeavyAtomFactory(RbtModelList modelList);
        virtual void VisitReceptorFlexData(RbtReceptorFlexData*);
        virtual void VisitLigandFlexData(RbtLigandFlexData*);
        virtual void VisitSolventFlexData(RbtSolventFlexData*);
        RbtAtomRList GetAtomList() const {return m_atomList;}
    private:
        RbtAtomRList m_atomList;
    };
    
  //Class type string
  static RbtString _CT;
  //Parameter names
  static RbtString _RMAX;
  static RbtString _QUADRATIC;//True = quadratic penalty function; false = linear
  
  RbtCavityGridSF(const RbtString& strName = "CAVITY");
  virtual ~RbtCavityGridSF();
  
 protected:
  virtual void SetupReceptor();
  virtual void SetupLigand();
  virtual void SetupSolvent();
  virtual void SetupScore();
  virtual RbtDouble RawScore() const;
  //DM 25 Oct 2000 - track changes to parameter values in local data members
  //ParameterUpdated is invoked by RbtParamHandler::SetParameter
  void ParameterUpdated(const RbtString& strName);
  
 private:
  RbtRealGridPtr m_spGrid;
  RbtDouble m_maxDist;//Max distance of any grid point from the cavity
  RbtAtomRList m_atomList;//combined list of all movable heavy atoms
  RbtDouble m_rMax;
  RbtBool m_bQuadratic;//synchronised with QUADRATIC named parameter
};

#endif //_RBTCAVITYGRIDSF_H_
