/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPMFDirSource.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
