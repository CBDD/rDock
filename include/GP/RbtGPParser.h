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

//GPParser Class. Parses a genome from a cartesian GP

#ifndef _RBT_GPPARSER_H_
#define _RBT_GPPARSER_H_

#include "RbtTypes.h"
#include "RbtContainers.h"
#include "RbtGPTypes.h"
#include "RbtCell.h"
#include "RbtCommand.h"
#include "RbtGPChromosome.h"

class RbtGPParser 
{
public:
    static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////
    RbtGPParser(RbtInt, RbtInt, RbtInt, RbtInt);
    RbtGPParser(const RbtGPParser&);
 
    ///////////////////
    // Destructor
    //////////////////
    virtual ~RbtGPParser();
  
	void Clear();
    ostream& Print(ostream&) const;
    friend ostream& operator<<(ostream& s, const RbtGPParser &p);
    RbtReturnTypeList Parse(RbtGPChromosomePtr,RbtReturnTypeList&);
    RbtString PrintParse(istream&,RbtGPChromosomePtr,RbtBool,RbtBool);
	RbtBool IsConstantInt(RbtInt i) { return (i == CINT);};
	RbtBool IsConstantFloat(RbtInt i) { return (i == CFLOAT);};

private:
    /////////////////////
    // Private data
    /////////////////////  
    RbtReturnType Eval(RbtGPChromosomePtr,RbtInt);
    RbtReturnType Parse1Output(RbtGPChromosomePtr,RbtInt);
    RbtString PrintEval(RbtGPChromosomePtr,RbtInt,RbtBool,RbtBool);
    RbtString PrintParse1Output(RbtGPChromosomePtr,RbtInt,RbtBool,RbtBool);
    //RbtGPChromosomePtr chrom;
    //RbtCellList cells;
    RbtCommandList commands;
    RbtInt nProgramInputs, nFunctionsInputs, nProgramOutputs;
	RbtInt ADD, SUB, MUL, DIV, IF, CINT, CFLOAT;
};

//Useful typedefs
typedef SmartPtr<RbtGPParser> RbtGPParserPtr;  //Smart pointer

#endif //_RbtGPParser_H_
