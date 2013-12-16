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

#include "RbtGPFitnessFunction.h"
#include "RbtGPGenome.h"
#include "RbtGPParser.h"
#include "RbtDebug.h"
#include <strstream>
#include <fstream>

RbtString RbtGPFitnessFunction::_CT("RbtGPFitnessFunction");

// Constructors
RbtGPFitnessFunction::RbtGPFitnessFunction() 
            : m_rand(Rbt::GetRbtRand()), objective(0.0), fitness(0.0)
{
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}
 
RbtGPFitnessFunction::RbtGPFitnessFunction(const RbtGPFitnessFunction& g) 
            : m_rand(Rbt::GetRbtRand()), objective(g.objective),
              fitness(g.fitness)
{
    _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtGPFitnessFunction& RbtGPFitnessFunction::operator=
                           (const RbtGPFitnessFunction& g)
{
    if (this != &g)
    {
        objective = g.objective;
        fitness = g.fitness;
    }
}


    ///////////////////
    // Destructor
    //////////////////
RbtGPFitnessFunction::~RbtGPFitnessFunction()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtDouble RbtGPFitnessFunction::GetObjective() const
{
    return objective;
}

RbtDouble RbtGPFitnessFunction::GetFitness() const
{
    return fitness;
}

void RbtGPFitnessFunction::SetFitness(RbtDouble f)
{
    fitness = f;
}

