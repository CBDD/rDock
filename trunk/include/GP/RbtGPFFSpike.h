/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtGPFFSpike.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

//GPFFSpike Class. GPFFSpike represents a fitness function
//valid for the CGP. It is used to train with the Spike set.

#ifndef _RBT_GPFFSPIKE_H_
#define _RBT_GPFFSPIKE_H_

#include "RbtGPTypes.h"
#include "RbtRand.h"
#include "RbtGPGenome.h"
#include "RbtGPFitnessFunction.h"
#include "RbtContext.h"
// #include <cmath>

class RbtGPFFSpike : public RbtGPFitnessFunction
{
public:
    static RbtString _CT;
    RbtGPFFSpike(RbtContextPtr c) : contextp(c)
	   {_RBTOBJECTCOUNTER_CONSTR_(_CT);};
    virtual ~RbtGPFFSpike()
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
    RbtGPFFSpike() ; // default constructor disabled

};

//Useful typedefs
typedef SmartPtr<RbtGPFFSpike> RbtGPFFSpikePtr;  //Smart pointer
typedef vector<RbtGPFFSpikePtr> RbtGPFFSpikeList;//Vector of smart pointers
typedef RbtGPFFSpikeList::iterator RbtGPFFSpikeListIter;
typedef RbtGPFFSpikeList::const_iterator RbtGPFFSpikeListConstIter;
#endif //_RbtGPFFSpike
