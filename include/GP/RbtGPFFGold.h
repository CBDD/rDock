/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtGPFFGold.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

//GPFFGold Class. GPFFGold represents a fitness function
//valid for the CGP. It is used to train with the Gold set.

#ifndef _RBT_GPFFGOLD_H_
#define _RBT_GPFFGOLD_H_

#include "RbtGPTypes.h"
#include "RbtRand.h"
#include "RbtGPGenome.h"
#include "RbtGPFitnessFunction.h"
// #include <cmath>

class RbtGPFFGold : public RbtGPFitnessFunction
{
public:
    static RbtString _CT;

    void ReadTables(istream&);
    RbtDouble CalculateFitness(RbtGPGenomePtr, RbtReturnTypeArray&,
                               RbtReturnTypeArray&, RbtDouble, RbtBool);
    RbtDouble CalculateFitness(RbtGPGenomePtr, RbtReturnTypeArray&,
                               RbtReturnTypeArray&, RbtBool);

};

//Useful typedefs
typedef SmartPtr<RbtGPFFGold> RbtGPFFGoldPtr;  //Smart pointer
typedef vector<RbtGPFFGoldPtr> RbtGPFFGoldList;//Vector of smart pointers
typedef RbtGPFFGoldList::iterator RbtGPFFGoldListIter;
typedef RbtGPFFGoldList::const_iterator RbtGPFFGoldListConstIter;
#endif //_RbtGPFFGold
