/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/GP/RbtGPGenome.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include "RbtGPGenome.h"
#include "RbtGPParser.h"
#include "RbtDebug.h"
#include <strstream>
#include <fstream>
using std::ifstream;
using std::ios_base;

RbtString RbtGPGenome::_CT("RbtGPGenome");
RbtInt RbtGPGenome::npi;
RbtInt RbtGPGenome::nfi;
RbtInt RbtGPGenome::nsfi;
RbtInt RbtGPGenome::no;
RbtInt RbtGPGenome::nn;
RbtInt RbtGPGenome::nf;
RbtInt RbtGPGenome::nr;
RbtInt RbtGPGenome::nc;
RbtInt RbtGPGenome::l;

// Constructors
RbtGPGenome::RbtGPGenome() : m_rand(Rbt::GetRbtRand())
{
    chrom = new RbtGPChromosome(npi, nfi, nn, no, nr, nc);
    fitness = 0.0;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}
 
RbtGPGenome::RbtGPGenome(const RbtGPGenome& g) : m_rand(Rbt::GetRbtRand())
{
    chrom = new RbtGPChromosome(*(g.chrom));
    fitness = g.fitness;
    _RBTOBJECTCOUNTER_COPYCONSTR_(_CT);
}

RbtGPGenome::RbtGPGenome(istream& in) : m_rand(Rbt::GetRbtRand())
{
    long int seed;
    in >> seed;
    m_rand.Seed(seed);
    // Get structure
    in >> npi >> nfi >> nsfi >> no >> nf >> nr >> nc >> l;
    nn  = nr * nc;    //number of nodes
    chrom = new RbtGPChromosome(npi, nfi, nn, no, nr, nc);
    in >> *chrom;
    fitness = 0.0;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtGPGenome::RbtGPGenome(RbtString str) : m_rand(Rbt::GetRbtRand())
{
    istrstream ist(str.c_str());
    // Get structure
    ist >> npi >> nfi >> nsfi >> no >> nf >> nr >> nc >> l;
    nn  = nr * nc;    //number of nodes
    chrom = new RbtGPChromosome(npi, nfi, nn, no, nr, nc);
    ist >> *chrom;
    fitness = 0.0;
}

RbtGPGenome& RbtGPGenome::operator=(const RbtGPGenome& g)
{
    if (this != &g)
    {
        *chrom = *(g.chrom);
        fitness = g.fitness;
    }
}


    ///////////////////
    // Destructor
    //////////////////
RbtGPGenome::~RbtGPGenome()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

void RbtGPGenome::SetStructure(RbtInt inpi, RbtInt infi, 
                               RbtInt insfi, RbtInt ino, 
                               RbtInt inf, RbtInt inr, 
                               RbtInt inc, RbtInt il)
{
    npi = inpi;   //number of program inputs
    nfi = infi;   //number of inputs per function
    nsfi = insfi; //number of inputs needed to calculate the SF
    no  = ino;    //number of outputs of the program. 
                  //Assumption: functions return one output
    nf  = inf;    //number of functions
    nr  = inr;    //number of rows
    nc  = inc;    //number of columns
    nn  = nr * nc;    //number of nodes
    l   = il;     //connectivity level, i.e., how many previous columns
                  //of cells may have their outputs connected to a
                  //node in the current column
}

void RbtGPGenome::Initialise()
{
    RbtInt s = chrom->size();
    RbtInt g = 0, max, min;
    for (RbtInt i = 0 ; i < nc ; i++)
    {
        for (RbtInt j = 0 ; j < nr ; j++)
        {
            for (RbtInt k = 0 ; k < nfi ; k++)
            {
                max = npi + i * nr;
                if (i < l)
                    (*chrom)[g++] = m_rand.GetRandomInt(max);
                else
                {
                    min = npi + (i - l) * nr;
                    (*chrom)[g++] = m_rand.GetRandomInt(max - min) + min;
                }
            }
            // node's function
            RbtInt f = m_rand.GetRandomInt(nf+1);
            if (f == nf)  // constant Int
            {
                RbtInt a, b;
                a = m_rand.GetRandomInt(200) - 100;
                b = m_rand.GetRandomInt(10) - 5;
                RbtDouble r = (a / 10.0) * pow(10.0,b);//m_rand.GetRandom01() * 20 - 10;
                chrom->SetConstant(r, g);
            }
            (*chrom)[g++] = f;
        }
    }
    //outputs of program
    min = npi + (nc - l) * nr;
    max = npi + nc * nr;
    for (RbtInt i = 0 ; i < no ; i++)
        (*chrom)[g++] = m_rand.GetRandomInt(max - min) + min;
    fitness = 0.0;
}

void RbtGPGenome::MutateGene(RbtInt i)
{
    RbtInt min, max, column;
    // is it part of the nodes?
    if (i < (nr * nc * (nfi + 1)))
    {
        // is it an input of a function?
        if (((i + 1) % (nfi + 1)) != 0)
        {
            column = i / (nr * (nfi + 1));
            max = npi + column * nr;
            if (column < l)
                (*chrom)[i] = m_rand.GetRandomInt(max);
            else
            {
                min = npi + (column - l) * nr;
                (*chrom)[i] = m_rand.GetRandomInt(max - min) + min;
            }
        }
        else // it is a function
        {
            RbtInt f = m_rand.GetRandomInt(nf+1); //m_rand.GetRandomInt(nf);
            if (f == nf)  // constant Int
            {
                RbtInt a, b;
                a = m_rand.GetRandomInt(200) - 100;
                b = m_rand.GetRandomInt(10) - 5;
                RbtDouble r = (a / 10.0) * pow(10.0,b);//m_rand.GetRandom01() * 20 - 10;
                chrom->SetConstant(r, i);
            }
            else if ((*chrom)[i] == nf)
            {
                chrom->ResetConstant(i);
            }
            (*chrom)[i] = f;
        }
    }
    else // it is an output
    {
        min = npi + (nc - l) * nr;
        max = npi + nc * nr;
        (*chrom)[i] = m_rand.GetRandomInt(max - min) + min;
    }
}

void RbtGPGenome::Mutate(RbtDouble mutRate)
{
    for (RbtInt i = 0; i < chrom->size() ; i++)
        if (m_rand.GetRandom01() < mutRate)
        {
            MutateGene(i);
        }
}


void RbtGPGenome::UniformCrossover(const RbtGPGenome& mom, 
                                   const RbtGPGenome& dad)
{
    RbtInt coin;
    for (RbtInt i = 0 ; i < mom.chrom->size() ; i++)
    {
        coin = m_rand.GetRandomInt(2);
        if (coin == 0)
            (*chrom)[i] = (*mom.chrom)[i];
        else
            (*chrom)[i] = (*dad.chrom)[i];
    }
}

void RbtGPGenome::Crossover(RbtGPGenome&){}

ostream& RbtGPGenome::Print(ostream& s) const
{
    RbtInt g = 0;
    // Print structure
    s << npi << " " << nfi << " " << nsfi << " " << no << " " <<
         nf  << " " << nr   << " " << nc << " " << l   << endl;
    s << *chrom << endl;
    RbtGPParser p(npi,nfi,0,no);
    ifstream nstr("descnames",ios_base::in);
    s << p.PrintParse(nstr, chrom, true, false);
    return s;
}

ostream& operator<<(ostream& s, const RbtGPGenome &g)
{
    return g.Print(s);
}

