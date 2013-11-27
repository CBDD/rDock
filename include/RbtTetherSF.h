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
