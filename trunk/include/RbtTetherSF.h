/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtTetherSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Tethered atoms distance restraint scoring function

#ifndef _RBTTETHERSF_H_
#define _RBTTETHERSF_H_

#include "RbtBaseInterSF.h"

class RbtTetherSF : public RbtBaseInterSF
{
    public:
    //Class type string
    static RbtString _CT;
    //Parameter names
    static RbtString _REFERENCE_FILE;
        
    RbtTetherSF(const RbtString& strName = "TETHER");
    virtual ~RbtTetherSF();
    
    protected:
    virtual void SetupReceptor();
    virtual void SetupLigand();
    virtual void SetupScore();
    virtual RbtDouble RawScore() const;
    //DM 25 Oct 2000 - track changes to parameter values in local data members
    //ParameterUpdated is invoked by RbtParamHandler::SetParameter
    void ParameterUpdated(const RbtString& strName);

    private:
	RbtIntList ReadTetherAtoms(RbtStringList &);
    RbtAtomList m_ligAtomList;
    RbtIntList m_tetherAtomList;
    RbtCoordList m_tetherCoords;
};

#endif //_RBTTETHERSF_H_
