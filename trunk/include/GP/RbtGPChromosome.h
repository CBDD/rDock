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

//GP chromosome. Keeps the list of integers that represent a chromosome and
//a list of cells that represent a program.

#ifndef _RBT_GPCHROMOSOME_H
#define _RBT_GPCHROMOSOME_H

#include "RbtGPTypes.h"
#include "RbtCell.h"
#include "RbtCommands.h"

class RbtGPChromosome 
{
public:
    static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////
    RbtGPChromosome(RbtInt, RbtInt, RbtInt, RbtInt, RbtInt, RbtInt);
    RbtGPChromosome(const RbtGPChromosome&);
 
    ///////////////////
    // Destructor
    //////////////////
    virtual ~RbtGPChromosome();
  
    void Clear();
    ostream& Print(ostream&) const;
    friend ostream& operator<<(ostream& s, const RbtGPChromosome &p);
    friend istream& operator>>(istream& s, RbtGPChromosome &p);
    const int& operator[] (int idx) const  {return chrom[idx];};
    int& operator[] (int idx) {return chrom[idx];};
    RbtGPChromosome& operator=(const RbtGPChromosome& c);
    RbtInt size() const {return chrom.size();};
    RbtCellPtr Cells(RbtInt idx) const
        {return (cells[idx]);}
    void SetConstant(RbtReturnType, RbtInt);
    void ResetConstant(RbtInt);
    RbtInt GetStartingCell() { return (chrom[chrom.size() - 1]); };
	RbtReturnType GetRCte(RbtInt ncell) { return cells[ncell]->GetResult();}; 
    RbtCommands GetCommand(RbtInt ncell) 
    {
        if (ncell < nProgramInputs)
            throw RbtError(_WHERE_, "This cell does not have a command");
        RbtCommands comm(chrom[(nFunctionsInputs + 1) * 
                               (ncell - nProgramInputs) + nFunctionsInputs]);
        return comm;
    }
    RbtInt GetArgument(RbtInt ncell, RbtInt narg)
    {
        if (ncell < nProgramInputs)
            throw RbtError(_WHERE_, "This cell does not have arguments");
        return (chrom[(nFunctionsInputs + 1) * 
                      (ncell - nProgramInputs) + narg]);
    }
    RbtBool IsProgramInput(RbtInt ncell)
    {
        return (ncell < nProgramInputs);
    }


private:
    /////////////////////
    // Private data
    /////////////////////  
    RbtIntList chrom;
    RbtCellList cells;
    RbtInt nProgramInputs, nFunctionsInputs, nProgramOutputs, nRows, 
           nColumns;
};

//Useful typedefs
typedef SmartPtr<RbtGPChromosome> RbtGPChromosomePtr;  //Smart pointer

#endif //_RBTGPCHROMOSOME_H_
