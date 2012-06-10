/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtGPParser.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
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
