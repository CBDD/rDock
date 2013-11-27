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

#include <sstream>
using std::ostringstream;
using std::ends;

#include "RbtCharmmDataSource.h"
#include "RbtCharmmTypesFileSource.h"
#include "RbtElementFileSource.h"
#include "RbtFileError.h"

////////////////////////////////////////
//Constructors/destructors

//DM 30 Apr 1999 - pass in masses.rtf file name as parameter (default=GetRbtFileName("data","masses.rtf"))
RbtCharmmDataSource::RbtCharmmDataSource(const RbtString& strMassesFile)
{
  //Need element data source to convert element names to atomic numbers
  RbtElementFileSourcePtr spElementData(new RbtElementFileSource(Rbt::GetRbtFileName("data","RbtElements.dat")));

  //Initialize m_atomTypes, m_atomicNumber and m_hybridState from masses.rtf
  RbtCharmmTypesFileSourcePtr spTypesFile(new RbtCharmmTypesFileSource(strMassesFile));
  RbtCharmmTypeList typeList = spTypesFile->GetTypeList();
  for (RbtCharmmTypeListIter iter = typeList.begin(); iter != typeList.end(); iter++) {
    m_atomTypes[iter->nAtomType] = iter->strAtomType;
    RbtElementData elData = spElementData->GetElementData(iter->element);
    m_atomicNumber[iter->strAtomType] = elData.atomicNo;
    //The hybridisation state has been encoded in the comments field in the Rbt version of masses.rtf
    m_hybridState[iter->strAtomType] = ConvertCommentStringToHybridState(iter->comment);
  }
					
  //Initialize m_implicitHydrogens
  m_implicitHydrogens["CH1E"] = 1;
  m_implicitHydrogens["CH2E"] = 2;
  m_implicitHydrogens["CH3E"] = 3;
  m_implicitHydrogens["C5RE"] = 1;
  m_implicitHydrogens["C6RE"] = 1;
  m_implicitHydrogens["SH1E"] = 1;
  m_implicitHydrogens["NH2E"] = 2; //Added DM 4 Jan 1999 - support older Charmm atom types in XPLOR

  //Initialize m_formalCharge
  //Note: these formal charges are not intended to sum to the correct overall charge on the molecule
  //It is only the sign of the charge which is correct, and indicates the atom can take part in
  //a charged interaction
  m_formalCharge["OC"] = -1; //e.g. in carboxylate, phosphate
  m_formalCharge["HC"] = +1; //e.g. in ammonium, guanidinium
  m_formalCharge["NC"] = +1; //e.g. in guanidinium

  _RBTOBJECTCOUNTER_CONSTR_("RbtCharmmDataSource");
}
//Default destructor
RbtCharmmDataSource::~RbtCharmmDataSource()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtCharmmDataSource");
}


////////////////////////////////////////
//Public methods
////////////////
RbtString RbtCharmmDataSource::AtomTypeString(RbtInt nAtomType) throw (RbtError)
{
 RbtIntStringMapIter iter = m_atomTypes.find(nAtomType);
  if (iter != m_atomTypes.end())
    return (*iter).second;
  else {
    ostringstream ostr;
    ostr << "CHARMm atom type #" << nAtomType << " not found in CharmmDataSource" << ends;
    RbtString strError(ostr.str());
    //delete ostr.str();
    throw RbtFileMissingParameter(_WHERE_,strError);
  }
}

RbtInt RbtCharmmDataSource::ImplicitHydrogens(const RbtString& strFFType) throw (RbtError)
{
 RbtStringIntMapIter iter = m_implicitHydrogens.find(strFFType);
  if (iter != m_implicitHydrogens.end())
    return (*iter).second;
  else
    return 0;
}

RbtInt RbtCharmmDataSource::AtomicNumber(const RbtString& strFFType) throw (RbtError)
{
 RbtStringIntMapIter iter = m_atomicNumber.find(strFFType);
  if (iter != m_atomicNumber.end())
    return (*iter).second;
  else
    throw RbtFileMissingParameter(_WHERE_,"CHARMm atom type " +strFFType+" not found in CharmmDataSource");
}

RbtInt RbtCharmmDataSource::FormalCharge(const RbtString& strFFType) throw (RbtError)
{
  RbtStringIntMapIter iter = m_formalCharge.find(strFFType);
  if (iter != m_formalCharge.end()) {
    return (*iter).second;
  }
  else
    return 0;
}

//DM 8 Dec 1998
RbtAtom::eHybridState RbtCharmmDataSource::HybridState(const RbtString& strFFType) throw (RbtError)
{
  RbtStringHybridStateMapIter iter = m_hybridState.find(strFFType);
  if (iter != m_hybridState.end())
    return (*iter).second;
  else
    throw RbtFileMissingParameter(_WHERE_,"CHARMm atom type " +strFFType+" not found in CharmmDataSource");
}

////////////////////////////////////////
//Private methods
////////////////////////////////////////

//DM 8 Dec 1998 Searches for Hybridisation state string in the masses.rtf comment field
//and returns the matching RbtAtom::eHybridState enum
//Valid strings are (RBT::SP), (RBT::SP2), (RBT::SP3), (RBT::AROM), (RBT::TRI)
//(brackets are important)
RbtAtom::eHybridState RbtCharmmDataSource::ConvertCommentStringToHybridState(const RbtString& strComment)
{
  if (strComment.find("(RBT::SP)") != RbtString::npos)
    return RbtAtom::SP;
  else if (strComment.find("(RBT::SP2)") != RbtString::npos)
    return RbtAtom::SP2;
  else if (strComment.find("(RBT::SP3)") != RbtString::npos)
    return RbtAtom::SP3;
  else if (strComment.find("(RBT::AROM)") != RbtString::npos)
    return RbtAtom::AROM;
  else if (strComment.find("(RBT::TRI)") != RbtString::npos)
    return RbtAtom::TRI;
  else
    return RbtAtom::UNDEFINED;
}
