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

#ifndef _RBTPMFDIRSOURCE_H_
#define _RBTPMFDIRSOURCE_H_

#include "RbtDirectorySource.h"
#include "RbtPMF.h"

class RbtPMFDirSource : public RbtDirectorySource
{
	public:
		RbtPMFDirSource(const RbtString&) throw (RbtError);
		virtual ~RbtPMFDirSource() {};

		virtual void 	ReadFiles(vector< vector<RbtPMFValue> > *aVect, vector<string> *aNameVect, vector<RbtPMFValue> *aSlopeVect);
		virtual void	ParseLines(vector<RbtString> anStrVect, vector<RbtPMFValue> *aValueVect);
};

#endif // _RBTPMFDIRSOURCE_H_
