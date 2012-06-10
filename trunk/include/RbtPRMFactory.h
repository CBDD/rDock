/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtPRMFactory.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Factory class for generating workspace objects from <receptor>.prm file
//Responsible for generation of:
//Receptor, Ligand, Solvent,
#ifndef RBTPRMFACTORY_H_
#define RBTPRMFACTORY_H_

#include "RbtModel.h"
#include "RbtError.h"
#include "RbtBaseMolecularFileSource.h"

class RbtParameterFileSource;
class RbtDockingSite;

class RbtPRMFactory {
    public:
        static const RbtString& _CT;
        //Receptor parameters
        static const RbtString& _REC_SECTION;
        static const RbtString& _REC_FILE;
        static const RbtString& _REC_TOPOL_FILE;
        static const RbtString& _REC_COORD_FILE;
        static const RbtString& _REC_NUM_COORD_FILES;
        static const RbtString& _REC_FLEX_DISTANCE;
        static const RbtString& _REC_DIHEDRAL_STEP;
        
        //Ligand parameters
        static const RbtString& _LIG_SECTION;
        //Ligand file is not read here, so no need for _LIG_FILE
        //No need to specify ligand flexibility parameter names explicitly
        //as we assume that all parameters in SECTION LIGAND
        //are flexibility params.
        
        //Solvent parameters
        static const RbtString& _SOLV_SECTION;
        static const RbtString& _SOLV_FILE;
        //No need to specify solvent flexibility parameter names explicitly
        //as we assume that all parameters in SECTION SOLVENT
        //are flexibility params (except for _SOLV_FILE)
  
        RbtPRMFactory(RbtParameterFileSource* pParamSource);
        RbtPRMFactory(RbtParameterFileSource* pParamSource,
                        RbtDockingSite* pDS);
        
        RbtInt GetTrace() const {return m_iTrace;}
        void SetTrace(RbtInt iTrace) {m_iTrace = iTrace;}
        RbtDockingSite* GetDockingSite() const {return m_pDS;}
        void SetDockingSite(RbtDockingSite* pDS) {m_pDS = pDS;}
       
        RbtModelPtr CreateReceptor() throw (RbtError);
        RbtModelPtr CreateLigand(RbtBaseMolecularFileSource* pSource) throw (RbtError);
        RbtModelList CreateSolvent() throw (RbtError);
        
    private:
        //Creates the appropriate source according to the file extension
        RbtMolecularFileSourcePtr CreateMolFileSource(
                const RbtString& fileName) throw (RbtError);
        void AttachReceptorFlexData(RbtModel* pReceptor);
        void AttachLigandFlexData(RbtModel* pLigand);
        void AttachSolventFlexData(RbtModel* pSolvent);
        RbtParameterFileSource* m_pParamSource;
        RbtDockingSite* m_pDS;
        RbtInt m_iTrace;
};
#endif /*RBTPRMFACTORY_H_*/
