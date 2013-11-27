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
