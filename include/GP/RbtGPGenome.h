/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/GP/RbtGPGenome.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

//GPGenome Class. GPGenome represents a genome in the CGP

#ifndef _RBT_GPGENOME_H_
#define _RBT_GPGENOME_H_

#include "RbtGPTypes.h"
#include "RbtRand.h"
#include "RbtGPChromosome.h"
// #include <cmath>

class RbtGPGenome 
{
public:
    static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////
    RbtGPGenome();
    RbtGPGenome(const RbtGPGenome&);
    RbtGPGenome(RbtString);
    RbtGPGenome(istream&);

    RbtGPGenome& operator=(const RbtGPGenome&);
 
    ///////////////////
    // Destructor
    //////////////////
    virtual ~RbtGPGenome();
    static void SetStructure(RbtInt, RbtInt, RbtInt, RbtInt,
                             RbtInt, RbtInt, RbtInt, RbtInt);
    static RbtInt GetNIP(){return npi;} ;
    static void SetNIP(RbtInt n){ npi = n;} ;
    static RbtInt GetNIF(){return nfi;} ;
    static RbtInt GetNN(){return nn;} ;
    static RbtInt GetNO(){return no;} ;
    static RbtInt GetNSFI(){return nsfi;} ;
    static void SetNSFI(RbtInt n){nsfi = n;} ;
    RbtGPChromosomePtr GetChrom() const {return (new RbtGPChromosome(*chrom));} ;
    void Initialise();
    void Mutate(RbtDouble);
	void SetFitness(RbtDouble f){fitness = f;};
	RbtDouble GetFitness()const{return fitness;};
    void UniformCrossover(const RbtGPGenome&, const RbtGPGenome&);
    void Crossover(RbtGPGenome&);
    friend ostream& operator<<(ostream& s, const RbtGPGenome &p);
    ostream& Print(ostream&) const;


protected:
  
private:
    void MutateGene(RbtInt);

    static RbtInt npi, nfi, nsfi, no, nn, nf, nr, nc, l;
        // defines structures and constraints of each genome
    RbtRand& m_rand;
    RbtDouble fitness;      // scaled score. For now same than score
    RbtGPChromosomePtr chrom;      // list of integers that represent a chrom
};

//Useful typedefs
typedef SmartPtr<RbtGPGenome> RbtGPGenomePtr;  //Smart pointer
typedef vector<RbtGPGenomePtr> RbtGPGenomeList;//Vector of smart pointers
typedef RbtGPGenomeList::iterator RbtGPGenomeListIter;
typedef RbtGPGenomeList::const_iterator RbtGPGenomeListConstIter;

#endif //_RbtGPGenome_H_
