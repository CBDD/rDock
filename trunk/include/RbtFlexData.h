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

//Abstract base class for all model flexibility types
//Provides Accept method for RbtFlexDataVisitors
//(subclasses must implement)
#ifndef _RBTFLEXDATA_
#define _RBTFLEXDATA_

#include "RbtFlexDataVisitor.h"
#include "RbtParamHandler.h"

class RbtModel;
class RbtDockingSite;

class RbtFlexData : public RbtParamHandler {
public:
	virtual ~RbtFlexData() {};
    RbtModel* GetModel() const {return m_pModel;}
    void SetModel(RbtModel* pModel) {m_pModel = pModel;}
    RbtDockingSite* GetDockingSite() const {return m_pDockSite;}
	virtual void Accept(RbtFlexDataVisitor&) = 0;
    
protected:
    RbtFlexData(RbtDockingSite* pDockSite)
            : m_pModel(NULL), m_pDockSite(pDockSite) {}

private:
	RbtModel* m_pModel;
    RbtDockingSite* m_pDockSite;
};

typedef SmartPtr<RbtFlexData> RbtFlexDataPtr;

#endif //_RBTFLEXDATA_
