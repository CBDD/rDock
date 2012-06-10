/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtCharmmDataSource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Class for providing Charmm force field data

#ifndef _RBTCHARMMDATASOURCE_H_
#define _RBTCHARMMDATASOURCE_H_

#include "RbtConfig.h"
#include "RbtAtom.h"

//Added by DM, 8 Dec 1998 for HybridState lookup method
typedef map<RbtString,RbtAtom::eHybridState> RbtStringHybridStateMap;
typedef RbtStringHybridStateMap::iterator RbtStringHybridStateMapIter;
typedef RbtStringHybridStateMap::const_iterator RbtStringHybridStateMapConstIter;

class RbtCharmmDataSource
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors

  //DM 30 Apr 1999 - pass in masses.rtf file name as parameter
  RbtCharmmDataSource(const RbtString& strMassesFile = Rbt::GetRbtFileName("data","masses.rtf"));

  ~RbtCharmmDataSource(); //Default destructor


  ////////////////////////////////////////
  //Public methods
  ////////////////
  RbtString AtomTypeString(RbtInt nAtomType) throw (RbtError);
  RbtInt ImplicitHydrogens(const RbtString& strFFType) throw (RbtError);
  RbtInt AtomicNumber(const RbtString& strFFType) throw (RbtError);
  RbtInt FormalCharge(const RbtString& strFFType) throw (RbtError);//DM 24 Mar 1999 - changed from double to int
  RbtAtom::eHybridState HybridState(const RbtString& strFFType) throw (RbtError);//DM 8 Dec 1998

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////


 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////

  RbtCharmmDataSource(const RbtCharmmDataSource&);//Copy constructor disabled by default

  RbtCharmmDataSource& operator=(const RbtCharmmDataSource&);//Copy assignment disabled by default

  //DM 8 Dec 1998 Searches for Hybridisation state string in the masses.rtf comment field
  //Valid strings are (RBT::SP), (RBT::SP2), (RBT::SP3), (RBT::AROM), (RBT::TRI)
  //(brackets are important)
  RbtAtom::eHybridState ConvertCommentStringToHybridState(const RbtString& strComment);

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtIntStringMap m_atomTypes;
  RbtStringIntMap m_implicitHydrogens;
  RbtStringIntMap m_atomicNumber;
  RbtStringIntMap m_formalCharge;
  RbtStringHybridStateMap m_hybridState;//DM 8 Dec 1998
};

//useful typedefs
typedef SmartPtr<RbtCharmmDataSource> RbtCharmmDataSourcePtr;//Smart pointer

#endif //_RBTCHARMMDATASOURCE_H_
