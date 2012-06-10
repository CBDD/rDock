/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtSolventFlexData.h#4 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Manages flexibility data for solvent
#ifndef _RBTSOLVENTFLEXDATA_
#define _RBTSOLVENTFLEXDATA_

#include "RbtLigandFlexData.h"

class RbtSolventFlexData : public RbtLigandFlexData {
public:
    //Probability of solvent being visible
    static const RbtString& _OCCUPANCY;
    //Nominal step size for mutating the occupancy variable
    static const RbtString& _OCCUPANCY_STEP;
    RbtSolventFlexData(RbtDockingSite* pDockSite);
	virtual void Accept(RbtFlexDataVisitor& v) {v.VisitSolventFlexData(this);}
     
private:
    RbtSolventFlexData();
};
#endif //_RBTSOLVENTFLEXDATA_
