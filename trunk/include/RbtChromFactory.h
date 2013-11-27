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
