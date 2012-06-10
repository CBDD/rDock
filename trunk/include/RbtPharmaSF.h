/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPharmaSF.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Pharmacophore distance restraint scoring function

#ifndef _RBTPHARMASF_H_
#define _RBTPHARMASF_H_

#include "RbtBaseInterSF.h"
#include "RbtConstraint.h"
#include "RbtBaseMolecularFileSink.h"

class RbtPharmaSF : public RbtBaseInterSF
{
    public:
    //Class type string
    static RbtString _CT;
    //Parameter names
    static RbtString _CONSTRAINTS_FILE;
    static RbtString _OPTIONAL_FILE;
    static RbtString _NOPT;
    static RbtString _WRITE_ERRORS;
        
    RbtPharmaSF(const RbtString& strName = "PHARMA");
    virtual ~RbtPharmaSF();
    //Override RbtBaseSF::ScoreMap to provide additional raw descriptors
    virtual void ScoreMap(RbtStringVariantMap& scoreMap) const;
    
    protected:
    virtual void SetupReceptor();
    virtual void SetupLigand();
    virtual void SetupScore();
    virtual RbtDouble RawScore() const;
    //DM 25 Oct 2000 - track changes to parameter values in local data members
    //ParameterUpdated is invoked by RbtParamHandler::SetParameter
    void ParameterUpdated(const RbtString& strName);

    private:
    RbtConstraintList m_constrList;
    RbtConstraintList m_optList;
    RbtInt m_nopt;
    RbtMolecularFileSinkPtr m_spErrorFile;
    RbtBool m_bWriteErrors;
    //Keep track of individual constraint scores for ScoreMap
    mutable RbtDoubleList m_conScores;//Mandatory constraint scores
    mutable RbtDoubleList m_optScores;//Optional constraint scores
};

#endif //_RBTPHARMASF_H_
