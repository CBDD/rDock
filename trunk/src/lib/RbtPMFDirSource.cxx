/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/lib/RbtPMFDirSource.cxx#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtPMFDirSource.h"
// #include <stdlib>
//#include "RbtFileError.h"

using std::cout;
using std::endl;

RbtPMFDirSource::RbtPMFDirSource(const RbtString& aDir) throw (RbtError): 
RbtDirectorySource(aDir)
{
#ifdef _DEBUG
	cout << _CT << " RbtPMFDirSource constructor" << endl;
#endif
}


void RbtPMFDirSource::ReadFiles(
				vector< vector<RbtPMFValue> > *aVect, 
				vector<string> *aNameVect,
				vector<RbtPMFValue> *aSlopeVect)
{
	cout<<"Reading "<< fNum <<" PMF files..."<<endl; 
	RbtString 		theFileName;		// will be the filename with full path
	RbtString 		theLine;			// one line from the file
	vector<string>	theStrData;			// data represented as strings in vector
	RbtPMFValue		theDummy = {0.0,0.0};	// dummy initial value for plateau 	
		
	while (fNum--) {
		theFileName += thePath+"/";
//		cout << "Does not work with Solaris" <<endl;
//		exit(1);
		RbtString theFileStr(fNameList[fNum]->d_name);				// copy from C-string 
//		RbtString theFileStr("junk");
		if ( theFileName.size()+theFileStr.size() > PATH_SIZE )	// check size (though should be ok)
		{	
			cout << _CT << theFileName << endl;
		   	cout << _CT << theFileStr << endl;
			throw RbtStringTooLong(_WHERE_,"");
		}
		else 
			theFileName += theFileStr;							// concatenate to get final size
		// get dot index in extension
		string::size_type theExtIdx = theFileName.find_last_of(".pmf");
		// this is for storing values from one file
		vector <RbtPMFValue>	theValues;
		// when file has no .pmf extension, skip. Also check for existence and perm
		if (string::npos != theExtIdx && !stat(theFileName.c_str(),&fStat) && S_ISREG(fStat.st_mode) ){ 
#ifdef _DEBUG
			cout << _CT << " Processing: " << theFileStr << endl;
#endif	// debug
			// get rid of .pmf and put it into the vector that will be used later to figure out types
			aNameVect->push_back(theFileStr.erase(theFileStr.find_last_of("."),4) );
			inFile.open(theFileName.c_str(),ifstream::in);
			// read file contents into vector
			while(inFile >> theLine)
				theStrData.push_back(theLine);
			inFile.close();
			ParseLines(theStrData, &theValues);
			aVect->push_back(theValues); 
			// search for the start of the 3.0 plateau backwards in pmf
			aSlopeVect->push_back(theDummy);	// create a dummy value
			for(int i=theValues.size();i>=0;i--) {
				if(theValues[i].density>=3.0) {
					aSlopeVect->back()= theValues[i];
					//aSlopeVect->back().distance = theValues[i].distance;
					//aSlopeVect->back().density  = theValues[i].density;
					break;
				}
			}
#ifdef GETJUNK
//#ifdef _DEBUG
			cout << _CT <<" Starting slope at: "<<aSlopeVect->back()<<" for ";
			cout <<theValues[i].density<<" ";
			cout <<theFileStr.substr(0,2)<<" "<<theFileStr.substr(2,2)<<endl; 
//#endif // _DEBUG
#endif //GETJUNK
			theStrData.erase(theStrData.begin(), theStrData.end()); // delete vector after parsing
		} 
#ifdef _DEBUG
		else {
			cout << _CT << " Skipping file: " << theFileName << endl;
		}
#endif	//_DEBUG
		
		theFileName.erase();
	}
	cout << " done."<<endl;
}

void RbtPMFDirSource::ParseLines(vector<RbtString> anStrVect, vector<RbtPMFValue> *aValueVect)
{
	RbtPMFValue				theActualValue;
	
	for (int i=0; i<anStrVect.size(); i++) {
		theActualValue.distance = (float) atof(anStrVect[i].data());
		i++;	// get next number
		theActualValue.density  = (float) atof(anStrVect[i].data());
		aValueVect->push_back(theActualValue);
	}
}

