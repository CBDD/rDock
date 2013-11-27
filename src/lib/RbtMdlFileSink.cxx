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

#include <iomanip>
using std::setw;
using std::setfill;

#include "RbtMdlFileSink.h"

////////////////////////////////////////
//Constructors/destructors
RbtMdlFileSink::RbtMdlFileSink(const RbtString& fileName, RbtModelPtr spModel) :
  RbtBaseMolecularFileSink(fileName,spModel),m_bFirstRender(true)
{
 //DM 27 Apr 1999 - default behaviour is for the first Render to overwrite any existing file
 //then subsequent Renders to append.
  SetAppend(false);
  //Open an Element data source
  m_spElementData = RbtElementFileSourcePtr(new RbtElementFileSource(Rbt::GetRbtFileName("data","RbtElements.dat")));
  _RBTOBJECTCOUNTER_CONSTR_("RbtMdlFileSink");
}

RbtMdlFileSink::~RbtMdlFileSink()
{
  _RBTOBJECTCOUNTER_DESTR_("RbtMdlFileSink");
}


////////////////////////////////////////
//Public methods
////////////////
//
////////////////////////////////////////
//Override public methods from RbtBaseFileSink
void RbtMdlFileSink::Render() throw (RbtError)
{
  try {
    RbtModelPtr spModel(GetModel());
    RbtAtomList modelAtomList(spModel->GetAtomList());
    RbtBondList modelBondList(spModel->GetBondList());
    RbtModelList solventList = GetSolvent();
    //Concatenate all solvent atoms and bonds into a single list
    //DM 7 June 2006 - only render the enabled solvent models
    RbtAtomList solventAtomList;
    RbtBondList solventBondList;
    for (RbtModelListConstIter iter = solventList.begin(); iter != solventList.end(); ++iter) {
          if ((*iter)->GetEnabled()) {
  	    RbtAtomList atomList = (*iter)->GetAtomList();
  	    RbtBondList bondList = (*iter)->GetBondList();
   	    std::copy(atomList.begin(),atomList.end(),std::back_inserter(solventAtomList));
  	    std::copy(bondList.begin(),bondList.end(),std::back_inserter(solventBondList));
          }
    }
    
    //Title lines
    RbtStringList titleList(spModel->GetTitleList());
    if (!titleList.empty())
      AddLine(titleList[0]);
    else
      AddLine(spModel->GetName());
    AddLine(RbtString("  rDOCK(R)          3D"));
    AddLine(Rbt::GetProduct()+"/"+Rbt::GetVersion()+"/"+Rbt::GetBuild());
    
    //Write number of atoms and bonds
    ostrstream ostr;
    ostr << setw(3) << modelAtomList.size() + solventAtomList.size()
	 << setw(3) << modelBondList.size() + solventBondList.size()
	 << setw(3) << 0 << setw(3) << 0 << setw(3) << 0
	 << setw(3) << 0 << setw(3) << 0 << setw(3) << 0
	 << setw(3) << 0 << setw(3) << 0 << setw(3) << 999 << " V2000"
	 << ends;
    AddLine(ostr.str());
    delete ostr.str();
    
    //DM 19 June 2006 - clear the map of logical atom IDs each time
    //we render a model
    m_atomIdMap.clear();
    RenderAtomList(modelAtomList);
    RenderAtomList(solventAtomList);
    RenderBondList(modelBondList);
    RenderBondList(solventBondList);
    RenderData(spModel->GetDataMap());
    
    AddLine("$$$$");
    
    Write();
    if (m_bFirstRender) {
      SetAppend(true);
      m_bFirstRender = false;
    }
  }
  catch (RbtError& error) {
    //Close();
    throw;//Rethrow the RbtError
  }
}
  
void RbtMdlFileSink::RenderAtomList(const RbtAtomList& atomList) {
	for (RbtAtomListConstIter aIter = atomList.begin(); aIter != atomList.end(); aIter++) {
	  RbtAtomPtr spAtom(*aIter);
	  //DM 19 June 2006. Check if this atom has been rendered previously.
	  if (m_atomIdMap.find(spAtom) == m_atomIdMap.end()) {
	    RbtUInt nextAtomId = m_atomIdMap.size()+1;
	    //cout << "RenderAtom " << spAtom->GetFullAtomName() << " (actual ID=" << spAtom->GetAtomId()
	    //	 << "); file ID=" << nextAtomId << endl;
	    m_atomIdMap.insert(make_pair(spAtom,nextAtomId));
	  }
	  else {
	    //Should never happen. Probably best to throw an error at this point.
	    throw RbtBadArgument(_WHERE_,"Duplicate atom rendered: " + spAtom->GetFullAtomName());
	  }
	  RbtElementData elData = m_spElementData->GetElementData(spAtom->GetAtomicNo());
	  RbtInt nFormalCharge = spAtom->GetFormalCharge();
	  if (nFormalCharge != 0)
		nFormalCharge = 4 - nFormalCharge;
	  ostrstream ostr;
	  ostr.precision(4);
	  ostr.setf(ios_base::fixed,ios_base::floatfield);
	  ostr.setf(ios_base::right,ios_base::adjustfield);
	  ostr << setw(10) << spAtom->GetX()
	   << setw(10) << spAtom->GetY()
	   << setw(10) << spAtom->GetZ();//X,Y,Z coord
	  ostr.setf(ios_base::left,ios_base::adjustfield);
	  ostr << setw(0) << " "
	   << setw(3) << elData.element.c_str();//Element name
	  ostr.setf(ios_base::right,ios_base::adjustfield);
	  ostr << setw(2) << 0 //mass difference
	   << setw(3) << nFormalCharge //charge
	   << setw(3) << 0 //atom stereo parity
	   << setw(3) << 0 //hydrogen count+1 (query CTABs only)
	   << setw(3) << 0 //stereo care box (query CTABs only)
	   << setw(3) << 0 //valence (0 = no marking)
	   << ends;//Mass diff, formal charge, stereo parity, num hydrogens, center
	  AddLine(ostr.str());
	  delete ostr.str();
	}
}

void RbtMdlFileSink::RenderBondList(const RbtBondList& bondList) {
	for (RbtBondListConstIter bIter = bondList.begin(); bIter != bondList.end(); bIter++) {
	  RbtBondPtr spBond(*bIter);
	  //DM 19 June 2006. Lookup the logical atom IDs for each atom in the bond
	  RbtAtomIdMap::const_iterator aIter1 = m_atomIdMap.find(spBond->GetAtom1Ptr());
	  RbtAtomIdMap::const_iterator aIter2 = m_atomIdMap.find(spBond->GetAtom2Ptr());
	  if ( (aIter1 != m_atomIdMap.end()) && (aIter2 != m_atomIdMap.end()) ) {
	    RbtUInt id1 = (*aIter1).second;
	    RbtUInt id2 = (*aIter2).second;
	    //cout << "RenderBond " << spBond->GetAtom1Ptr()->GetFullAtomName()
	    //	 << spBond->GetAtom2Ptr()->GetFullAtomName()
	    //	 << "; file ID1=" << id1
	    //	 << "; file ID2=" << id2 << endl;
	    ostrstream ostr;
	    ostr.setf(ios_base::right,ios_base::adjustfield);
	    ostr << setw(3) << id1
		 << setw(3) << id2
		 << setw(3) << spBond->GetFormalBondOrder()
		 << setw(3) << 0
		 << setw(3) << 0
		 << setw(3) << 0
		 << ends;//Atom1, Atom2, bond order, stereo designator, unused, topology code
	    AddLine(ostr.str());
	    delete ostr.str();
	  }
	  else {
	    //Should never happen. Probably best to throw an error at this point.
	    throw RbtBadArgument(_WHERE_,"Error rendering bond, logical atom IDs not found"); 
	  }
	}
}

void RbtMdlFileSink::RenderData(const RbtStringVariantMap& dataMap) {
	if (!dataMap.empty()) {
		AddLine("M  END");//End of molecule marker
		for (RbtStringVariantMapConstIter iter = dataMap.begin(); iter != dataMap.end(); iter++) {
			AddLine(">  <"+(*iter).first+">");//Field name
			RbtStringList sl = (*iter).second.StringList();
			for (RbtStringListConstIter slIter = sl.begin(); slIter != sl.end(); ++slIter) {
				AddLine(*slIter);//Field values
			}
			AddLine("");//Blank line denotes end of field values
		}
	}
}

