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

#include "RbtGPPopulation.h"
#include "RbtDebug.h"
#include <algorithm>
#include <float.h>

      // Static data members for RbtGPPopulation
RbtString RbtGPPopulation::_CT("RbtGPPopulation"); 


      // Constructor
      
RbtGPPopulation::RbtGPPopulation(RbtInt s, RbtInt nr, 
                              RbtGPFitnessFunctionPtr f,
                              RbtReturnTypeArray& it,
                              RbtReturnTypeArray& sft) 
                            : m_rand(Rbt::GetRbtRand())
{
    popsize = s;
    nrepl = nr;
    ff = f;
    ittrain = it;
    sfttrain = sft;
    cout << "const: " << ittrain.size() <<endl;
    c = 2.0;        // default value for the sigma truncation multiplier
    pop = RbtGPGenomeList(popsize);
    for (RbtGPGenomeListIter iter = pop.begin() ; iter != pop.end() ; iter++)
        *iter = RbtGPGenomePtr(new RbtGPGenome());
    newpop = RbtGPGenomeList (nrepl);
    for (RbtGPGenomeListIter iter = newpop.begin() ; 
                             iter != newpop.end() ; iter++)
        *iter = RbtGPGenomePtr(new RbtGPGenome());
    psum = new RbtDouble [popsize];
    _RBTOBJECTCOUNTER_CONSTR_(_CT);

}

      // Copy constructor
      // This is copying the pointers. Do I need to change it to 
      // copying the actual genomes?
RbtGPPopulation::RbtGPPopulation(const RbtGPPopulation & p)
                                : m_rand(Rbt::GetRbtRand())
{
    popsize = p.popsize;
    nrepl = p.nrepl;
    ff = p.ff;
    ittrain = p.ittrain;
    sfttrain = p.sfttrain;
    c = 2.0;        // default value for the sigma truncation multiplier
    pop = RbtGPGenomeList(popsize);
    for (RbtInt i = 0 ; i < p.popsize ; i++)
        pop[i] = p.pop[i];
    newpop = RbtGPGenomeList (p.newpop.size());
    for (RbtInt i = 0 ; i < p.newpop.size() ; i++)
        newpop[i] = p.newpop[i];
    psum = new RbtDouble [popsize];
    _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

      // Destructor
      
RbtGPPopulation::~RbtGPPopulation()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

RbtInt RbtGPPopulation::GetSize()
{
    return popsize;
}

RbtInt RbtGPPopulation::GetNrepl()
{
    return nrepl;
}

       // Initialises the population by calling the initializer for each genome
void RbtGPPopulation::Initialise(RbtDouble hitlimit,RbtBool function)
{
    RbtDouble bestFitness = -FLT_MAX, fit;
    for (RbtInt i = 0 ; i < popsize ; i++)
    {
        pop[i]->Initialise(); 
        fit = ff->CalculateFitness(pop[i],ittrain,sfttrain,function);
        cout << "init: " << ittrain.size() << endl;
        if (fit > bestFitness)
        {
            bestFitness = fit;
            bestInd = i;
        }
    }
}

// calculates the scoring function. This is normally calculated as soon as a new
// indidivual is created. But we need to make a call to this function when the 
// slope is changed, to make sure all the individuals get the new value for 
// their scoring function before continuing with the GA
void RbtGPPopulation::Eval(RbtDouble hitlimit, RbtBool function)
{
    for (RbtInt i = 0 ; i < popsize ; i++)
        ff->CalculateFitness(pop[i],ittrain,sfttrain,function);
}

// Size 2 probabilistic tournament selection, where the winner of the
// tournament is accepted with a probability of tp. Otherwise the
// loser is accepted.
RbtGPGenomePtr RbtGPPopulation::TSelect(RbtDouble tp) const
{
    RbtInt ind1, ind2;
    ind1 = m_rand.GetRandomInt(popsize);
    do
    {
        ind2 = m_rand.GetRandomInt(popsize);
    } while (ind1 == ind2);
    RbtDouble r = m_rand.GetRandom01();
    if (pop[ind1]->GetFitness() > pop[ind2]->GetFitness())
        // winner is ind1
        if (r < tp)
            return (pop[ind1]);
        else
            return (pop[ind2]);
    else
        // winner is ind2
        if (r < tp)
            return (pop[ind2]);
        else
            return (pop[ind1]);
}



 
RbtGPGenomePtr RbtGPPopulation::Select(RbtString selector) const
{
    return TSelect(0.7);
}

void RbtGPPopulation::GAstep(RbtString selector, RbtDouble pcross,
                             RbtDouble pmut, 
                             RbtDouble mean, RbtDouble variance,
                             RbtDouble hitlimit, RbtBool function) 
                                throw (RbtError)
{
    RbtGPGenomePtr mother, father;
    for (RbtInt i = 0 ; i < nrepl-1 ; i++)
    {
        mother = Select(selector);
        father = Select(selector);
        RbtInt j = 0;
        while (father == mother)
        {
            father = Select(selector);
            if (j > 100)  
              // only if all or almost all individuals of the population are 
              // copies of the same one
                  throw RbtError(_WHERE_, 
                          "Population failure - not enough diversity");
            j++;
        }
        if (m_rand.GetRandom01() < pcross)
            newpop[i]->UniformCrossover(*mother, *father);
        newpop[i]->Mutate(pmut);
    }
    *(newpop[nrepl-1]) = *(pop[bestInd]);

      // calculate the objective values and 
      // sort newpop  newpop = RbtGPGenomeList (p.newpop.size());
    RbtDouble bestFitness = -FLT_MAX, fit;
    for (RbtInt i = 0 ; i < newpop.size() ; i++)
    {
        *(pop[i]) = *(newpop[i]);
        fit = ff->CalculateFitness(pop[i],ittrain,sfttrain,function);
        if (fit > bestFitness)
        {
            bestFitness = fit;
            bestInd = i;
        }
    }
}

void RbtGPPopulation::EPstep(RbtString selector, RbtDouble pcross,
                             RbtDouble pmut, 
                             RbtDouble mean, RbtDouble variance,
                             RbtDouble hitlimit, RbtBool function) 
                                throw (RbtError)
{
    *(newpop[0]) = *(pop[bestInd]);
    for (RbtInt i = 1 ; i < popsize ; i++)
    {
        *(newpop[i]) =  *(pop[bestInd]);
//        cout << "new ind\n";
        newpop[i]->Mutate(pmut);
    }
      // calculate the objective values and 
      // sort newpop  newpop = RbtGPGenomeList (p.newpop.size());
    RbtDouble bestFitness = -FLT_MAX, fit;
    for (RbtInt i = 0 ; i < newpop.size() ; i++)
    {
        *(pop[i]) = *(newpop[i]);
        fit = ff->CalculateFitness(pop[i],ittrain,sfttrain,function);
        if (fit >= bestFitness)
        {
            bestFitness = fit;
            bestInd = i;
        }
    }
}


// the worst individual from the population are replaced by the best from p
RbtGPGenomePtr RbtGPPopulation::Best() const
{
    return (pop[bestInd]);
}

ostream& RbtGPPopulation::Print(ostream& s) const
{
    for (RbtInt i = 0 ; i < pop.size() ; i++)
    {
        pop[i]->Print(s);
        s << pop[i]->GetFitness() << endl;
    }
    return s;
}

ostream& operator<<(ostream&s, const RbtGPPopulation& p)
{
    return p.Print(s);
}

// calculates the scaled objectives 
// of the current population. To calculate the scaled objectives, we 
// use sigma truncating scaling as explained in Golberg's book,
// page 124
void RbtGPPopulation::ScaleFitness()
{ 
    total = 0;
      // calculate the new Fitnesses and the average 
    for (int i = 0 ; i < popsize ; i++)
    {
        total += pop[i]->GetFitness();
    }
    ave = total / popsize;  
      // calculate variance
    RbtDouble ss = 0.0;
    for (int i = 0 ; i < popsize ; i++)
    {
        ss += (pop[i]->GetFitness() - ave) * (pop[i]->GetFitness() - ave);
    }
      // standard deviation
    stdev = sqrt(ss / (popsize - 1));
      // calculate scaled values using sigma truncation
      // Goldberg page 124
    for (int i = 0 ; i < popsize ; i++)
    {
        RbtDouble f = pop[i]->GetFitness() - (ave - c * stdev);
        if (f < 0.0)
            f = 0.0;
        pop[i]->SetFitness(f);
    }
}

RbtGPGenomePtr RbtGPPopulation::RwSelect() const
{
    RbtDouble cutoff;
    RbtInt i, upper, lower;

    cutoff = m_rand.GetRandom01();
    lower = 0;
    upper = popsize - 1;
    while (upper >= lower)
    {
        i = lower + (upper - lower) / 2;
        if (psum[i] > cutoff)
            upper = i - 1;
        else
            lower = i + 1;
    }
      // make sure lower is a number between 0 and size - 1
    lower = (lower > (popsize - 1)) ? (popsize - 1) : lower;
    lower = (lower < 0) ? 0 : lower;
    return (pop[lower]);
}


void RbtGPPopulation::Swap( RbtGPPopulationPtr p, RbtInt nmigr)
{
    for (RbtInt i = 0 ; i < nmigr ; i++)
    {
       RbtGPGenomeListIter where = find(pop.begin(), pop.end(), p->pop[i]);
       //if (!where)
       if (where == pop.end())
          // copy it only if it is not already there
            this->pop[popsize - nmigr + i] = p->pop[i];
    }
}

RbtGPGenomePtr RbtGPPopulation::RkSelect() const
{
    RbtDouble bias = 1.9;
    RbtDouble index;
    RbtInt max = popsize; // index will be between 0 and popsize - 1
    index = max * (bias - sqrt(bias * bias - 4.0 
                * (bias - 1) * m_rand.GetRandom01()))
                / 2.0 / (bias - 1);
    if (index < 0)
    {
        cout << index << " " << max << endl;
        index = 0.0;
    }
    //DM 25/04/05 - avoid compiler warning by explicitly converting double to int
    RbtInt intIndex = (int)floor(index);
    return (pop[intIndex]);
}

/*
RbtBool RbtGPPopulation::Gen_eq(RbtGPGenome* x, RbtGPGenome* y)
{
    return true; //x->GenEqual(*x, *y);
}

void RbtGPPopulation::QSort(RbtGPGenomeList &p)
{
     sort(p.begin(), p.end(), GenCmp());
}

void RbtGPPopulation::QSort()
{
    sort(pop.begin(), pop.end(), GenCmp());
}

void RbtGPPopulation::MergePops()
{
    RbtGPGenomeList mpop(popsize + newpop.size(), pop[0]);
#ifdef _DEBUG
    cout << "pop\n";
    for (RbtInt k1 = 0 ; k1 < popsize ; k1++)
        cout << *(pop[k1]) << "  " << pop[k1]->GetFitness() << "  " 
             << pop[k1]->GetFitness() << endl;
    cout << "newpop\n";
    for (RbtGPGenomeListIter p = newpop.begin() ; p != newpop.end() ; p++)
        cout << (*p) << "  " << (*p)->GetFitness() 
             << "  " << (*p)->GetFitness() << endl;
#endif
    merge(pop.begin(), pop.end(), newpop.begin(), newpop.end(), 
          mpop.begin(), GenCmp());
    RbtGPGenomeListIter end = unique(mpop.begin(), mpop.end(), Gen_eq); 
            // uses operator== , Gen_eq);
    mpop.erase(end, mpop.end());
#ifdef _DEBUG
    cout << "mpop\n";
    for (RbtInt k1 = 0 ; k1 < mpop.size() ; k1++)
        cout << mpop[k1] << "  " << mpop[k1]->GetFitness() << "  " 
             << mpop[k1]->GetFitness() << endl;
#endif
    pop.clear();
    end = (mpop.size() > popsize) ? (mpop.begin() + popsize) : mpop.end();
    copy(mpop.begin(), end, back_inserter(pop));
#ifdef _DEBUG
    cout << "copied pop\n";
    for (RbtInt k1 = 0 ; k1 < popsize ; k1++)
        cout << *(pop[k1]) << "  " << pop[k1]->GetFitness() << "  " 
             << pop[k1]->GetFitness() << endl;
#endif
}*/
      
