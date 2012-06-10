/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtFlexData.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
