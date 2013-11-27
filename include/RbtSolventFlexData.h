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
