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
#include <iomanip>
using std::setw;
#include <stdio.h>
#include "RbtPsfFileSink.h"

////////////////////////////////////////
//Constructors/destructors

RbtPsfFileSink::RbtPsfFileSink(const RbtString& fileName, RbtModelPtr spModel) :
  RbtBaseMolecularFileSink(fileName,spModel) //Call base class constructor
{
  SetMultiConf(false);//Default to allowing multiple conformers/ models in crd files
  _RBTOBJECTCOUNTER_CONSTR_("RbtPsfFileSink");
}

RbtPsfFileSink::~RbtPsfFileSink()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtPsfFileSink");
}


////////////////////////////////////////
//Public methods
////////////////
//
////////////////////////////////////////
//Override public methods from RbtBaseFileSink
void RbtPsfFileSink::Render() throw (RbtError)
{
  //String constants in PSF files
  const RbtString strPsfKey("PSF");
  const RbtString strTitleKey("!NTITLE");
  const RbtString strAtomKey("!NATOM");
  const RbtString strBondKey("!NBOND:");

  //Reset() has the intelligence to only reset the counters when required
  //so we can simply call it each time here
  Reset();

  try {
    RbtModelPtr spModel(GetModel());
    
    //1. Write the titles
    AddLine(strPsfKey);
    AddLine("");
    RbtStringList titleList(spModel->GetTitleList());
    ostringstream ostr;
    ostr << setw(8) << titleList.size() << setw(0) << " " << strTitleKey << ends;
    AddLine(ostr.str());
    //delete ostr.str();
    for (RbtStringListConstIter iter=titleList.begin(); iter != titleList.end(); iter++) {
			//Add an initial star if the title does not already commence with a star
			if ((*iter).find("*") != 0)
	  		AddLine("*"+(*iter));
			else
	  		AddLine(*iter);
		}
		AddLine("");
    
    //2. Write number of atoms
    ostringstream ostr2;
    ostr2 << setw(8) << spModel->GetNumAtoms() << setw(0) << " " << strAtomKey << ends;
    AddLine(ostr2.str());
    //delete ostr2.str();
    
    //3. Write the atom list
    RbtAtomList atomList(spModel->GetAtomList());
    char* line = new char[128];
    for (RbtAtomListConstIter aIter = atomList.begin(); aIter != atomList.end(); aIter++) {
      RbtAtomPtr spAtom(*aIter);

      m_nAtomId++;
      sprintf(line,"%8d %-4.4s %-4.4s %-4.4s %-4.4s %-4.4s%15.3f%10.3f%12.1f",
	      spAtom->GetAtomId(),
	      spAtom->GetSegmentName().c_str(),
	      spAtom->GetSubunitId().c_str(),
	      spAtom->GetSubunitName().c_str(),
	      spAtom->GetAtomName().c_str(),
	      spAtom->GetFFType().c_str(),
	      spAtom->GetPartialCharge(),
	      spAtom->GetAtomicMass(),
	      0.0);
      AddLine(line);
      //Render the atom to a string stream
      //ostringstream ostr;
      //ostr.precision(5);
      //ostr.setf(ios_base::fixed,ios_base::floatfield);
      //ostr.setf(ios_base::right,ios_base::adjustfield);
      //ostr.fill('#');
			//ostr << setw(8) << spAtom->GetAtomId() << " ";
      //ostr.setf(ios_base::left,ios_base::adjustfield);
      //ostr << setw(4) << spAtom->GetSegmentName() << " ";
      //ostr << setw(4) << spAtom->GetSubunitId() << " ";
			//ostr << setw(4) << spAtom->GetSubunitName() << " ";
	   	//ostr << setw(4) << spAtom->GetAtomName() << " ";
	   	//ostr << setw(4) << spAtom->GetFFType() << " ";
      //ostr.setf(ios_base::right,ios_base::adjustfield);
      //ostr << setw(14) << spAtom->GetPartialCharge();
      //ostr << setw(10) << spAtom->GetAtomicMass();
      //ostr << setw(12) << 0.0 << ends;
      //AddLine(ostr.str());
      //delete ostr.str();
		}
		delete [] line;
    //4. Write number of bonds
    ostringstream ostr3;
    ostr3 << setw(8) << spModel->GetNumBonds() << setw(0) << " " << strBondKey << ends;
		AddLine("");
    AddLine(ostr3.str());
    //delete ostr3.str();
    
    //5. Write the bond list
    RbtBondList bondList(spModel->GetBondList());
    RbtInt nDiv = bondList.size() / 4;
    RbtInt nRem = bondList.size() % 4;
    RbtBondListConstIter bIter = bondList.begin();
    for (RbtInt i = 0; i < nDiv; i++) {
    	ostringstream ostr;
    	ostr << setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
						<< setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
			bIter++;
    	ostr << setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
						<< setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
			bIter++;
    	ostr << setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
						<< setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
			bIter++;
    	ostr << setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
						<< setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
			bIter++;
			ostr << ends;
			AddLine(ostr.str());
			//delete ostr.str();
		}
		//Remaining bonds on final incomplete line
		if (nRem>0) {
			ostringstream ostr;
    	for (RbtInt i = 0; i < nRem; i++,bIter++) {
    		ostr << setw(8) << (*bIter)->GetAtom1Ptr()->GetAtomId()
							<< setw(8) << (*bIter)->GetAtom2Ptr()->GetAtomId();
			}
			ostr << ends;
			AddLine(ostr.str());
			//delete ostr.str();
		}
    
		Write();//Commit the cache to the file
	}
  catch (RbtError& error) {
    throw;//Rethrow the RbtError
  }
}
