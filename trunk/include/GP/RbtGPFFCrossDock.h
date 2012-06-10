/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtGPFFCrossDock.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

//GPFFCrossDock Class. GPFFCrossDock represents a fitness function
//valid for the CGP. It is used to train with the CrossDock set.

#ifndef _RBT_GPFFCROSSDOCK_H_
#define _RBT_GPFFCROSSDOCK_H_

#include "RbtGPTypes.h"
#include "RbtRand.h"
#include "RbtGPGenome.h"
#include "RbtGPFitnessFunction.h"
// #include <cmath>

class RbtGPFFCrossDock : public RbtGPFitnessFunction
{
public:
    static RbtString _CT;

    void ReadTables(istream&,RbtReturnTypeArray&,RbtReturnTypeArray&);
    RbtDouble CalculateFitness(RbtGPGenomePtr,RbtReturnTypeArray&,
			                   RbtReturnTypeArray&,RbtDouble,RbtBool);
    RbtDouble CalculateFitness(RbtGPGenomePtr,RbtReturnTypeArray&,
			                   RbtReturnTypeArray&, RbtBool);
private:
	void CreateRandomCtes(RbtInt);
	RbtDoubleList ctes;

};

//Useful typedefs
typedef SmartPtr<RbtGPFFCrossDock> RbtGPFFCrossDockPtr;  //Smart pointer
typedef vector<RbtGPFFCrossDockPtr> RbtGPFFCrossDockList;//Vector of smart pointers
typedef RbtGPFFCrossDockList::iterator RbtGPFFCrossDockListIter;
typedef RbtGPFFCrossDockList::const_iterator RbtGPFFCrossDockListConstIter;
#endif //_RbtGPFFCrossDock
