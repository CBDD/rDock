/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtGPFitnessFunction.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

//GPFitnessFunction Class. GPFitnessFunction represents a fitness function
//valid for the CGP

#ifndef _RBT_GPFITNESSFUNCTION_H_
#define _RBT_GPFITNESSFUNCTION_H_

#include "RbtGPTypes.h"
#include "RbtRand.h"
#include "RbtGPGenome.h"
#include "RbtContext.h"
// #include <cmath>

class RbtGPFitnessFunction 
{
public:
    static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////
    RbtGPFitnessFunction();
    RbtGPFitnessFunction(const RbtGPFitnessFunction&);

    RbtGPFitnessFunction& operator=(const RbtGPFitnessFunction&);
 
    ///////////////////
    // Destructor
    //////////////////
    virtual ~RbtGPFitnessFunction();
    virtual void ReadTables(istream&,RbtReturnTypeArray&,RbtReturnTypeArray&)=0;
    RbtDouble GetObjective() const;
    RbtDouble GetFitness() const;
    void  SetFitness(RbtDouble);
    virtual RbtDouble CalculateFitness(RbtGPGenomePtr,RbtReturnTypeArray&,
                                     RbtReturnTypeArray&, RbtDouble,RbtBool)=0;
    virtual RbtDouble CalculateFitness(RbtGPGenomePtr,RbtReturnTypeArray&,
                                     RbtReturnTypeArray&, RbtBool)=0;


protected:
  
    RbtRand& m_rand;
    RbtDouble objective;    // raw score
    RbtDouble fitness;      // scaled score
    RbtReturnTypeArray inputTable, SFTable;
private:
};

//Useful typedefs
typedef SmartPtr<RbtGPFitnessFunction> RbtGPFitnessFunctionPtr;  //Smart pointer
typedef vector<RbtGPFitnessFunctionPtr> RbtGPFitnessFunctionList;//Vector of smart pointers
typedef RbtGPFitnessFunctionList::iterator RbtGPFitnessFunctionListIter;
typedef RbtGPFitnessFunctionList::const_iterator RbtGPFitnessFunctionListConstIter;
#endif //_RbtGPFitnessFunction
