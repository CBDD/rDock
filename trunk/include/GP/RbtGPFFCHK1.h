/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtGPFFCHK1.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

//GPFFCHK1 Class. GPFFCHK1 represents a fitness function
//valid for the CGP. It is used to train with the CHK1 set.

#ifndef _RBT_GPFFCHK1_H_
#define _RBT_GPFFCHK1_H_

#include "RbtGPTypes.h"
#include "RbtRand.h"
#include "RbtGPGenome.h"
#include "RbtGPFitnessFunction.h"
#include "RbtContext.h"
// #include <cmath>

class RbtGPFFCHK1 : public RbtGPFitnessFunction
{
public:
    static RbtString _CT;
    RbtGPFFCHK1(RbtContextPtr c) : contextp(c)
	   {_RBTOBJECTCOUNTER_CONSTR_(_CT);};
    virtual ~RbtGPFFCHK1()
	   {_RBTOBJECTCOUNTER_DESTR_(_CT);};

    void ReadTables(istream&,RbtReturnTypeArray&,RbtReturnTypeArray&);
    RbtDouble CalculateFitness(RbtGPGenomePtr,RbtReturnTypeArray&,
                               RbtReturnTypeArray&,RbtDouble,RbtBool);
    RbtDouble CalculateFitness(RbtGPGenomePtr,RbtReturnTypeArray&,
                               RbtReturnTypeArray&, RbtBool);
private:
    void CreateRandomCtes(RbtInt);
    RbtDoubleList ctes;
    RbtCellContextPtr contextp;
    RbtGPFFCHK1() ; // default constructor disabled
		RbtReturnTypeArray inputTestTable, inputTrainingTable;
		RbtReturnTypeArray sfTestTable, sfTrainingTable;

};

//Useful typedefs
typedef SmartPtr<RbtGPFFCHK1> RbtGPFFCHK1Ptr;  //Smart pointer
typedef vector<RbtGPFFCHK1Ptr> RbtGPFFCHK1List;//Vector of smart pointers
typedef RbtGPFFCHK1List::iterator RbtGPFFCHK1ListIter;
typedef RbtGPFFCHK1List::const_iterator RbtGPFFCHK1ListConstIter;
#endif //_RbtGPFFCHK1
