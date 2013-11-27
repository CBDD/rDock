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

#include "RbtNmrRestraintFileSource.h"
#include "RbtFileError.h"

////////////////////////////////////////
//Constructors/destructors
RbtNmrRestraintFileSource::RbtNmrRestraintFileSource(const RbtString& fileName) :
  RbtBaseFileSource(fileName)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtNmrRestraintFileSource");
}

//Destructor
RbtNmrRestraintFileSource::~RbtNmrRestraintFileSource()
{
  ClearRestraintCache();
  _RBTOBJECTCOUNTER_DESTR_("RbtNmrRestraintFileSource");
}


////////////////////////////////////////
//Public methods
////////////////

//Number of NOE restraints in file
RbtUInt RbtNmrRestraintFileSource::GetNumNoeRestraints()
{
  Parse();
  return m_noeRestraintList.size();
}

//List of NOE restraints
RbtNoeRestraintNamesList RbtNmrRestraintFileSource::GetNoeRestraintList()
{
  Parse();
  return m_noeRestraintList;
}

//Number of STD restraints in file
RbtUInt RbtNmrRestraintFileSource::GetNumStdRestraints()
{
  Parse();
  return m_stdRestraintList.size();
}

//List of STD restraints
RbtStdRestraintNamesList RbtNmrRestraintFileSource::GetStdRestraintList()
{
  Parse();
  return m_stdRestraintList;
}  


////////////////////////////////////////
//Private methods
/////////////////

//Pure virtual in RbtBaseFileSource - needs to be defined here
void RbtNmrRestraintFileSource::Parse() throw (RbtError)
{
  const RbtString strDelimiter(",");//DM 23 Sept 1999

  //Only parse if we haven't already done so
  if (!m_bParsedOK) {
    ClearRestraintCache();//Clear current cache
    Read();//Read the file

    try {
      for (RbtFileRecListIter fileIter = m_lineRecs.begin(); fileIter != m_lineRecs.end(); fileIter++) {
	//Ignore blank lines and comment lines
	if ( ((*fileIter).length() == 0) || ((*fileIter).at(0) == '#') ) {
	  continue;
	}
	//Read the restraint
	RbtString strAtomNames1;
	RbtString strAtomNames2;
	RbtDouble maxDist(0.0);
	istrstream istr((*fileIter).c_str());
	istr >> strAtomNames1 >> strAtomNames2 >> maxDist;
	//Check if we read all the fields OK
	if (!istr)
	  throw RbtFileParseError(_WHERE_,"Missing field(s) at line "+(*fileIter)+" in "+GetFileName());

	//DM 10 Dec 2002 - support for STD restraints. The first atom name string is replaced by the
	//keyword STD (or std) in the file. The restraint is applied to the second atom list and implies
	//that the specified atoms should be within maxDist of *any* part of the receptor surface
	if (strAtomNames1 == "std" || strAtomNames1 == "STD") {
	  //Add new restraint to the list
	  RbtStdRestraintNames restraint;
	  restraint.from.type = NoeRestraintType(strAtomNames2);
	  //Check for bad syntax (restraint type=UNDEFINED)
	  if (restraint.from.type==Rbt::NOE_UNDEFINED)
	    throw RbtFileParseError(_WHERE_,"Unmatched or misplaced brackets at line "+(*fileIter)+" in "+GetFileName());
	  RbtStringList fromList = Rbt::ConvertDelimitedStringToList(strAtomNames2,strDelimiter);
	  restraint.from.names = fromList;
	  restraint.maxDist = maxDist;
	  m_stdRestraintList.push_back(restraint);
	}

	//NOE restraints
	else {
	  //Add new restraint to the list
	  RbtNoeRestraintNames restraint;
	  restraint.from.type = NoeRestraintType(strAtomNames1);
	  restraint.to.type = NoeRestraintType(strAtomNames2);
	  //Check for bad syntax (restraint type=UNDEFINED)
	  if ( (restraint.from.type==Rbt::NOE_UNDEFINED) || (restraint.to.type==Rbt::NOE_UNDEFINED))
	    throw RbtFileParseError(_WHERE_,"Unmatched or misplaced brackets at line "+(*fileIter)+" in "+GetFileName());
	  RbtStringList fromList = Rbt::ConvertDelimitedStringToList(strAtomNames1,strDelimiter);
	  restraint.from.names = fromList;
	  RbtStringList toList = Rbt::ConvertDelimitedStringToList(strAtomNames2,strDelimiter);
	  restraint.to.names = toList;
	  restraint.maxDist = maxDist;
	  m_noeRestraintList.push_back(restraint);
	}
      }
      //////////////////////////////////////////////////////////
      //If we get this far everything is OK
      m_bParsedOK = true;
    }
    
    catch (RbtError& error) {
      ClearRestraintCache();
      throw;//Rethrow the RbtError
    }
  }
}


void RbtNmrRestraintFileSource::ClearRestraintCache()
{
  m_noeRestraintList.clear();
  m_stdRestraintList.clear();
}

//Returns NOE restraint type and modifies the atom name string accordingly
//Returns UNDEFINED if the atom name string has bad syntax
Rbt::eNoeType RbtNmrRestraintFileSource::NoeRestraintType(RbtString& strAtomNames)
{
  Rbt::eNoeType restraintType(Rbt::NOE_UNDEFINED);

  //Restraint type: NOE_MEAN
  //Syntax: (at1,at2,at3...)
  //Check if atom names are bracketed correctly, if so change type to MEAN, and remove the brackets
  RbtString::size_type ob = strAtomNames.find("(");
  RbtString::size_type cb = strAtomNames.find(")");
  RbtString::size_type iLast = strAtomNames.length()-1;
  if ( (ob==0) && (cb==iLast)) {
    strAtomNames.erase(iLast,1);
    strAtomNames.erase(0,1);
    restraintType = Rbt::NOE_MEAN;
  }
  else if ( (ob==RbtString::npos) && (cb==RbtString::npos)) {
    //Restraint type: NOE_AND
    //Syntax: [at1,at2,at3...]
    //Check if atom names are bracketed correctly, if so change type to AND, and remove the brackets
    RbtString::size_type osqb = strAtomNames.find("[");
    RbtString::size_type csqb = strAtomNames.find("]");
    if ( (osqb==0) && (csqb==iLast)) {
      strAtomNames.erase(iLast,1);
      strAtomNames.erase(0,1);
      restraintType = Rbt::NOE_AND;
    }
    //Restraint type: NOE_OR
    //Syntax: at1,at2,at3...
    //Check for no brackets present
    else if ( (osqb==RbtString::npos) && (csqb==RbtString::npos)) {
      restraintType = Rbt::NOE_OR;
    }
  }
  return restraintType;
}
