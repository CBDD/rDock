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
