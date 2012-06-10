/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtGPFFHSP90.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

//GPFFHSP90 Class. GPFFHSP90 represents a fitness function
//valid for the CGP. It is used to train with the HSP90 set.

#ifndef _RBT_GPFFHSP90_H_
#define _RBT_GPFFHSP90_H_

#include "RbtGPTypes.h"
#include "RbtRand.h"
#include "RbtGPGenome.h"
#include "RbtGPFitnessFunction.h"
#include "RbtContext.h"

class RbtGPFFHSP90 : public RbtGPFitnessFunction
{
public:
    static RbtString _CT;
    RbtGPFFHSP90(RbtContextPtr c) : contextp(c)
	   {_RBTOBJECTCOUNTER_CONSTR_(_CT);};
    virtual ~RbtGPFFHSP90()
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
    RbtGPFFHSP90() ; // default constructor disabled
		RbtReturnTypeArray inputTestTable, inputTrainingTable;
		RbtReturnTypeArray sfTestTable, sfTrainingTable;

};

//Useful typedefs
typedef SmartPtr<RbtGPFFHSP90> RbtGPFFHSP90Ptr;  //Smart pointer
typedef vector<RbtGPFFHSP90Ptr> RbtGPFFHSP90List;//Vector of smart pointers
typedef RbtGPFFHSP90List::iterator RbtGPFFHSP90ListIter;
typedef RbtGPFFHSP90List::const_iterator RbtGPFFHSP90ListConstIter;
#endif //_RbtGPFFHSP90
