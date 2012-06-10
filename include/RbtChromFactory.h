/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtChromFactory.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Chromosome factory class
//Creates an RbtChrom aggregate by visiting a number of
//RbtFlexData subclasses.
//It is the users responsibility to delete the
//chromosome returned by GetChrom() when no longer required.
//No smart pointers are used within RbtChromFactory.
#ifndef _RBTCHROMFACTORY_H_
#define _RBTCHROMFACTORY_H_

#include "RbtFlexDataVisitor.h"
#include "RbtModelMutator.h"

class RbtChromElement;

class RbtChromFactory : public RbtFlexDataVisitor {
	public:
        RbtChromFactory();
		virtual void VisitReceptorFlexData(RbtReceptorFlexData*);
		virtual void VisitLigandFlexData(RbtLigandFlexData*);
        virtual void VisitSolventFlexData(RbtSolventFlexData*);

        RbtChromElement* GetChrom() const;
        //Temporary solution whilst we replace the ModelMutator class
        RbtModelMutatorPtr GetModelMutator() const {return m_spMutator;}
        
    private:
        RbtChromElement* m_pChrom;
        RbtModelMutatorPtr m_spMutator;
};
#endif //_RBTCHROMFACTORY_H_
