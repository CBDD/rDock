/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtFlexDataVisitor.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
