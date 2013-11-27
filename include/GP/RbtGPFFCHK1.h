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
