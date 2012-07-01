/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

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

