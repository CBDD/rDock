/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/GP/RbtGPFitnessFunction.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
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

