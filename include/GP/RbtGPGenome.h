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
