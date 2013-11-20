/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtCrdFileSink.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include <iomanip>
using std::setw;
using std::setfill;

#include "RbtCrdFileSink.h"

////////////////////////////////////////
//Constructors/destructors
//RbtCrdFileSink::RbtCrdFileSink(const char* fileName, RbtModelPtr spModel) :
//  RbtBaseMolecularFileSink(fileName,spModel) //Call base class constructor
//{
//  SetMultiConf(true);//Default to allowing multiple conformers/ models in crd files
//  _RBTOBJECTCOUNTER_CONSTR_("RbtCrdFileSink");
//}

RbtCrdFileSink::RbtCrdFileSink(const RbtString& fileName, RbtModelPtr spModel) :
  RbtBaseMolecularFileSink(fileName,spModel) //Call base class constructor
{
  SetMultiConf(true);//Default to allowing multiple conformers/ models in crd files
  _RBTOBJECTCOUNTER_CONSTR_("RbtCrdFileSink");
}

RbtCrdFileSink::~RbtCrdFileSink()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtCrdFileSink");
}


////////////////////////////////////////
//Public methods
////////////////
//
////////////////////////////////////////
//Override public methods from RbtBaseFileSink
void RbtCrdFileSink::Render() throw (RbtError)
{
  //Reset() has the intelligence to only reset the counters when required
  //so we can simply call it each time here
  Reset();

  try {
    RbtModelPtr spModel(GetModel());
    
    //Only write titles and number of atoms if this is the first conformer
    if (m_nAtomId == 0) {
      //1. Write the titles
      RbtStringList titleList(spModel->GetTitleList());
      //DM 8 June 1999 - limit number of title lines output to a maximum of 1
      RbtStringListConstIter iter;
      RbtInt nLines;
      for (iter=titleList.begin(),nLines=0; (iter != titleList.end()) && (nLines<1); iter++,nLines++) {
	//Add an initial star if the title does not already commence with a star
	if ((*iter).find("*") != 0)
	  AddLine("*"+(*iter));
	else
	  AddLine(*iter);
      }
      //Title delimiter record (* on its own)
      AddLine("*");

      //2. Write number of atoms
      ostrstream ostr;
      ostr << setw(5) << spModel->GetNumAtoms() << ends;
      AddLine(ostr.str());
      delete ostr.str();//DM 24 Mar 1999 - fix memory leak
      //Remember the line number containing the number of atoms
      //as we'll need to update the total num atoms after each Render when we are in multiconf mode
      //m_numAtomsLineRec = titleList.size()+1;
      m_numAtomsLineRec = nLines+1;
    }
    
    //In multiconf mode, we need to keep track of the current subunit and segment IDs so we can
    //detect when they change
    RbtString strLastSubunitId("UNLIKELY_TO_MATCH_THIS_BY_CHANCE");
    RbtString strLastSegmentName("UNLIKELY_TO_MATCH_THIS_BY_CHANCE");

    //3. Write the atom list
    RbtAtomList atomList(spModel->GetAtomList());
    for (RbtAtomListConstIter aIter = atomList.begin(); aIter != atomList.end(); aIter++) {
      RbtAtomPtr spAtom(*aIter);

      m_nAtomId++;
      //New subunit detected
      if (spAtom->GetSubunitId() != strLastSubunitId) {
	m_nSubunitId++;
	strLastSubunitId = spAtom->GetSubunitId();
      }
      //New segment detected
      if (spAtom->GetSegmentName() != strLastSegmentName) {
	m_nSegmentId++;
	strLastSegmentName = spAtom->GetSegmentName();
      }

      //Render the atom to a string stream
      ostrstream ostr;
      ostr.precision(5);
      ostr.setf(ios_base::fixed,ios_base::floatfield);
      ostr.setf(ios_base::right,ios_base::adjustfield);
      //In multiconf mode, use the incremental atom ID from the sink
      if (GetMultiConf())
	ostr << setw(5) << m_nAtomId;
      //Else use the atom ID from the atom itself
      else
	ostr << setw(5) << spAtom->GetAtomId();
      ostr << setw(5) << m_nSubunitId //This subunit ID is unique in the file
	   << setw(0) << " ";
      ostr.setf(ios_base::left,ios_base::adjustfield);
      ostr << setw(4) << spAtom->GetSubunitName().c_str()
	   << setw(0) << " "
	   << setw(4) << spAtom->GetAtomName().c_str();
      ostr.setf(ios_base::right,ios_base::adjustfield);
      ostr << setw(10) << spAtom->GetX()
	   << setw(10) << spAtom->GetY()
	   << setw(10) << spAtom->GetZ()
	   << setw(0) << " ";

      //DM 21 Apr 1999 - if we are not using the model segment names, use the incremental segment ID from the sink
      //DM 11 Feb 1999 - output leading zeros for segID, so they sort correctly when read as text
      if (!m_bUseModelSegmentNames) {
	ostr << setw(3) << setfill('0') << m_nSegmentId << setfill(' ') << " ";
	ostr.setf(ios_base::left,ios_base::adjustfield);
      }
      //Else use the segment name from the atom itself
      else {
	ostr.setf(ios_base::left,ios_base::adjustfield);
	ostr << setw(4) << spAtom->GetSegmentName().c_str();
      }
      
      ostr << setw(0) << " "
	   << setw(4) << spAtom->GetSubunitId().c_str(); //This subunit ID is unique in the segment
      ostr.setf(ios_base::right,ios_base::adjustfield);
      ostr << setw(10) << 0.0
	   << ends;
      AddLine(ostr.str());
      delete ostr.str();//DM 24 Mar 1999 - fix memory leak
    }

    //In multiconf mode, we need to update the total number of atoms in the cache
    //then write the file without clearing the cache
    if (GetMultiConf()) {
      ostrstream ostr;
      ostr << setw(5) << m_nAtomId << ends;
      ReplaceLine(ostr.str(),m_numAtomsLineRec);
      delete ostr.str();//DM 24 Mar 1999 - fix memory leak
      Write(false);
    }
    else
      Write();//Commit the cache to the file
  }
  catch (RbtError& error) {
    //Close();
    throw;//Rethrow the RbtError
  }
}
