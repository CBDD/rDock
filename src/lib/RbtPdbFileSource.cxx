/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtPdbFileSource.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

	#include "RbtPdbFileSource.h"
	
	//Constructors
	RbtPdbFileSource::RbtPdbFileSource(const char* fileName) :
	  RbtBaseMolecularFileSource(fileName,"PDB_FILE_SOURCE") //Call base class constructor
	{
	  _RBTOBJECTCOUNTER_CONSTR_("RbtPdbFileSource");
	}
	
	RbtPdbFileSource::RbtPdbFileSource(const RbtString fileName) :
	  RbtBaseMolecularFileSource(fileName,"PDB_FILE_SOURCE") //Call base class constructor
	{
	  _RBTOBJECTCOUNTER_CONSTR_("RbtPdbFileSource");
	}
	
	
	//Default destructor
	RbtPdbFileSource::~RbtPdbFileSource()
	{
	  _RBTOBJECTCOUNTER_DESTR_("RbtPdbFileSource");
	}
	
	void RbtPdbFileSource::Parse() throw (RbtError)
	{
	  //Expected string constants in PDB files
	  const RbtString strTitleKey("REMARK ");
	  const RbtString strAtomKey("ATOM ");
	  const RbtString strHetAtmKey("HETATM ");
	  //const RbtString strEndKey("END");
	 
	  //Only parse if we haven't already done so
	  if (!m_bParsedOK) {
	    ClearMolCache();//Clear current cache
	    Read();//Read the file
	
	    try {
	      for (RbtFileRecListIter fileIter = m_lineRecs.begin(); fileIter != m_lineRecs.end(); fileIter++) {
		//Ignore blank lines
		if ( (*fileIter).length() == 0 ) {
		  continue;
		}
		//Check for Title record
		else if ( (*fileIter).find(strTitleKey) == 0) {
		  RbtString strTitle = *fileIter;
		  strTitle.erase(0,strTitleKey.length());
		  m_titleList.push_back(strTitle);
		}
		//Check for Atom or HetAtm records
		else if ( ((*fileIter).find(strAtomKey) == 0) || ((*fileIter).find(strHetAtmKey) == 0) ) {
		  RbtString strDummy;
		  RbtInt nAtomId(1); //atom number
		  RbtString strAtomName; //atom name
	      RbtString strSegmentName; //segment name (chain identifier)
		  RbtString strSubunitName; //subunit(residue) name
		  RbtString strSubunitId; //subunit(residue) ID
		  RbtCoord coord; //X,Y,Z coords
		  //DM 15 June 2006 - change default occupancy to 1
		  //to avoid disabling solvent if pdb file does not contain
		  //occupancy field
		  RbtDouble occupancy(1.0);
		  RbtDouble tempFactor(0.0);
		  
		  RbtFileRec::size_type length = (*fileIter).size();
		  if (length > 10) {
		  	istringstream istr((*fileIter).substr(6,5).c_str());
		  	istr >> nAtomId;
			//cout << "Atom ID      =" << nAtomId << endl;
		  }
		  if (length > 15) {
		  	istringstream istr((*fileIter).substr(12,4).c_str());
		  	istr >> strAtomName;
			//cout << "Atom name    =" << strAtomName << endl;
		  }
		  if (length > 19) {
		    istringstream istr((*fileIter).substr(17,3).c_str());
		    istr >> strSubunitName;
			//cout << "Subunit name =" << strSubunitName << endl;
		  }
		  if (length > 21) {
		    istringstream istr((*fileIter).substr(21,1).c_str());
		    istr >> strSegmentName;
			//cout << "Segment name =" << strSegmentName << endl;
		  }
		  if (length > 25) {
		    istringstream istr((*fileIter).substr(22,4).c_str());
		    istr >> strSubunitId;
			//cout << "Subunit ID   =" << strSubunitId << endl;
		  }
		  if (length > 53) {
		    istringstream istr((*fileIter).substr(30,24).c_str());
		    istr >> coord.x
		  	   	 >> coord.y
		  		 >> coord.z;
			//cout << "coord        =" << coord << endl;
		  }
		  if (length > 59) {
		    istringstream istr((*fileIter).substr(54,6).c_str());
		    istr >> occupancy;
		    //cout << "occupancy    =" << occupancy << endl;
		  }
		  if (length > 65) {
		    istringstream istr((*fileIter).substr(60,6).c_str());
		    istr >> tempFactor;
		    //cout << "tempFactor   =" << tempFactor << endl;
		  }

		  //Construct a new atom (constructor only accepts the 2D params)
		  RbtAtomPtr spAtom(new RbtAtom(nAtomId,
						0,//Atomic number undefined
						strAtomName,
						strSubunitId,
						strSubunitName,
						strSegmentName
						)
				    );
	
	
		  //Now set the 3-D params we have
		  spAtom->SetCoords(coord);
          
          //Store the occupancy and temperature factor values
          //in User1 and User2 for use by solvent flexibility code
          //User1 is also used by scoring functions so these
          //values must be used immediately after creating the model
	      spAtom->SetUser1Value(occupancy);
          spAtom->SetUser2Value(tempFactor);
          
		  m_atomList.push_back(spAtom);
		  m_segmentMap[strSegmentName]++;//increment atom count in segment map
		}
	
		//////////////////////////////////////////////////////////
		//If we get this far everything is OK
		m_bParsedOK = true;
	      }
	    }
	
	    catch (RbtError& error) {
	      ClearMolCache();//Clear the molecular cache so we don't return incomplete atom and bond lists
	      throw;//Rethrow the RbtError
	    }
	  }
	}
