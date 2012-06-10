/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtTriposAtomType.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Atom typing class for Tripos 5.2 forcefield

#ifndef _RBTTRIPOSATOMTYPE_H_
#define _RBTTRIPOSATOMTYPE_H_

#include "RbtConfig.h"

class RbtAtom;

class RbtTriposAtomType
{
 public:
  ///////////////////////////////////////////////
  //Enums
  ///////////////////////////////////////////////

  //Atom types
  enum eType
  {
    UNDEFINED = 0,
    Al,
    Br,
    C_cat,
    C_1,C_1_H1,
    C_2,C_2_H1,C_2_H2,
    C_3,C_3_H1,C_3_H2,C_3_H3,
    C_ar,C_ar_H1,
    Ca,
    Cl,
    Du,
    F,
    H,H_P,
    I,
    K,
    Li,
    LP,
    N_1, N_2, N_3, N_4, N_am, N_ar, N_pl3,
    Na,
    O_2, O_3, O_co2,
    P_3,
    S_2, S_3, S_o, S_o2,
    Si,
    MAXTYPES //KEEP AS LAST TYPE: used to size the atom name string list
  };

  struct info {
    RbtString name;
    RbtInt atomicNo;
    int hybrid;
    info() : atomicNo(0), hybrid(0) {};
    info(const RbtString& n, RbtInt a, RbtInt h) : name(n), atomicNo(a), hybrid(h) {};
  };
      
  ////////////////////////////////////////
  //Constructors/destructors
  
  RbtTriposAtomType();
  virtual ~RbtTriposAtomType();
  
  ////////////////////////////////////////
  //Public methods
  ////////////////
  
  //Returns Tripos atom type for given atom
  //If useExtendedTypes is true, pseudo types for extended carbons and polar Hs are used
  eType operator() (RbtAtom* pAtom, RbtBool useExtendedTypes=false) const;
  //Converts Tripos type to string
  RbtString Type2Str(eType) const;
  // Get hybridisation from Tripos type
  int Type2Hybrid(eType) const;
  // Get atomic number from Tripos type
  RbtInt Type2AtomicNo(eType) const;
  //Converts string to Tripos type
  eType Str2Type(const RbtString&) const;
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  
  
 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  void SetupTypeInfo();
  
 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////
  
  
 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  vector<info> m_typeInfo;
  
};

//Useful typedefs
typedef vector<RbtTriposAtomType::info> RbtTriposInfoList;
typedef RbtTriposInfoList::iterator RbtTriposInfoListIter;
typedef RbtTriposInfoList::const_iterator RbtTriposInfoListConstIter;

typedef vector<RbtTriposAtomType::eType> RbtTriposAtomTypeList;
typedef RbtTriposAtomTypeList::iterator RbtTriposAtomTypeListIter;
typedef RbtTriposAtomTypeList::const_iterator RbtTriposAtomTypeListConstIter;

#endif //_RBTTRIPOSATOMTYPE_H_
