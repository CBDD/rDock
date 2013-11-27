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

//Abstract interface for all visitors of RbtFlexData subclasses
#ifndef _RBTFLEXDATAVISITOR_H_
#define _RBTFLEXDATAVISITOR_H_

class RbtReceptorFlexData;
class RbtLigandFlexData;
class RbtSolventFlexData;

class RbtFlexDataVisitor {
	public:
        virtual ~RbtFlexDataVisitor() {};
		virtual void VisitReceptorFlexData(RbtReceptorFlexData*) = 0;
		virtual void VisitLigandFlexData(RbtLigandFlexData*) = 0;
        virtual void VisitSolventFlexData(RbtSolventFlexData*) = 0;
};
#endif //_RBTFLEXDATAVISITOR_H_
