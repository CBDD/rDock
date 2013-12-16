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

#include "RbtGPParser.h"
#include "RbtCommand.h"
#include "RbtCell.h"
#include <strstream>
#include <fstream>
using std::ostrstream;

RbtString RbtGPParser::_CT("RbtGPParser");
RbtInt RbtCommand::ntabs = 0;
RbtBool RbtCommand::inside = false;

// Constructors
RbtGPParser::RbtGPParser(RbtInt npi, RbtInt nfi, RbtInt nn, RbtInt no)
{
    nProgramInputs = npi;
    nFunctionsInputs = nfi;
    nProgramOutputs = no;
    commands.push_back(new AddCommand());
    commands.push_back(new SubCommand());
    commands.push_back(new MulCommand());
    commands.push_back(new DivCommand());
    commands.push_back(new IfCommand());
	IF = 4;
    commands.push_back(new LogCommand());
    commands.push_back(new ExpCommand());
}

RbtGPParser::RbtGPParser(const RbtGPParser& p)
{
    nProgramInputs = p.nProgramInputs;
    nFunctionsInputs = p.nFunctionsInputs;
    nProgramOutputs = p.nProgramOutputs;
    commands.push_back(new AddCommand());
    commands.push_back(new SubCommand());
    commands.push_back(new MulCommand());
    commands.push_back(new DivCommand());
    commands.push_back(new IfCommand());
	IF = 4;
    commands.push_back(new LogCommand());
    commands.push_back(new ExpCommand());
}

///////////////////
// Destructor
//////////////////
RbtGPParser::~RbtGPParser(){}

RbtReturnType RbtGPParser::Eval(RbtGPChromosomePtr chrom, RbtInt n)
{
    if (!(chrom->Cells(n)->Evaluated()))
    {
        RbtReturnType o = Parse1Output(chrom, n);
        chrom->Cells(n)->SetResult(o);
        return o;
    }
    else
        return chrom->Cells(n)->GetResult();
}

RbtReturnType RbtGPParser::Parse1Output(RbtGPChromosomePtr chrom, RbtInt output)
{
    RbtInt startCell = (nFunctionsInputs + 1) * (output - nProgramInputs);
    RbtInt ncomm = (*chrom)[startCell + nFunctionsInputs];
    RbtInt nargs = commands[ncomm]->GetNArgs();
    RbtReturnType args[nargs];
    // I had this two loops as one before, but because is a
    // recursive function, they MUST be separated, otherwise the command
    // can get the arguments from another cell 
    for (RbtInt i = 0 ; i < nargs ; i++)
        args[i] = Eval(chrom, (*chrom)[startCell + i]);
    for (RbtInt i = 0 ; i < nargs ; i++)
        commands[ncomm]->SetArg(i, args[i]);
    RbtReturnType f2 = commands[ncomm]->Execute();
    return (commands[ncomm]->Execute());
}

RbtReturnTypeList RbtGPParser::Parse(RbtGPChromosomePtr chrom, 
                                   RbtReturnTypeList& inputs)
{
    // Give the values of the inputs to the corresponding cells
    for (RbtInt i = 0 ; i < nProgramInputs ; i++)
        chrom->Cells(i)->SetResult(*(inputs[i]));
    // for each output, find the node where its output is taken
    // from, and parse it
    RbtInt size = chrom->size();
    RbtReturnTypeList outputs;
    for (RbtInt i = 0 ; i < nProgramOutputs ; i++)
    {
        RbtReturnType o = Eval(chrom, 
                                       (*chrom)[size - nProgramOutputs + i]);
        outputs.push_back(new RbtReturnType(o));
    }
    return outputs;

}

RbtString RbtGPParser::PrintEval(RbtGPChromosomePtr chrom, RbtInt n, 
                               RbtBool numbers, RbtBool ins)
{
    if (!(chrom->Cells(n)->Named()))
    {
        RbtString o(PrintParse1Output(chrom, n,numbers,ins));
        ostrstream nm;
        if (numbers)
        {
            nm << "(" ;
            nm << n << " ";
            nm << o << ") " << ends;
        }
        else
            nm << o << ends;
        chrom->Cells(n)->SetName(nm.str());
        return (nm.str());
    }
    else
        return chrom->Cells(n)->GetName();
}

RbtString RbtGPParser::PrintParse1Output(RbtGPChromosomePtr chrom, RbtInt output,
                                       RbtBool numbers, RbtBool ins)
{
    RbtInt startCell = (nFunctionsInputs + 1) * (output - nProgramInputs);
    RbtInt ncomm = (*chrom)[startCell + nFunctionsInputs];
    if (ncomm == IF)
        RbtCommand::ntabs++;
    RbtInt nargs =  commands[ncomm]->GetNArgs();
    RbtStringList fs(nargs);
    for (RbtInt i = 0 ; i < commands[ncomm]->GetNArgs() ; i++)
        if ((ncomm == IF) && (i == 0))
        {
 //           RbtCommand::inside = true;
            fs[i] = PrintEval(chrom, (*chrom)[startCell + i], numbers, true);
//            RbtCommand::inside = false;
        }
        else
            fs[i] = PrintEval(chrom, (*chrom)[startCell + i], numbers,ins);
    for (RbtInt i = 0 ; i < commands[ncomm]->GetNArgs() ; i++)
        commands[ncomm]->SetNameArg(i, fs[i]);
    RbtString c;
    RbtCommand::inside = ins;
    if (ncomm == IF)
    {
        c = commands[ncomm]->GetName();
        RbtCommand::ntabs--;
    }
    else
        c = commands[ncomm]->GetName();
    return (c);
}

RbtString RbtGPParser::PrintParse(istream& istr, RbtGPChromosomePtr chrom, 
                                RbtBool numbers, RbtBool ins)
{
    // Give the names of the inputs to the corresponding cells
    RbtString name;
    RbtInt node;
    while (istr >> node)
    {
        istr.get();
        getline(istr, name);
        chrom->Cells(node)->SetName(name + " ");
    }
    // for each output, find the node where its output is taken
    // from, and parse it
    RbtInt size = chrom->size();
    RbtReturnTypeList outputs;
    RbtString f;
    for (RbtInt i = 0 ; i < nProgramOutputs ; i++)
        f += PrintEval(chrom, (*chrom)[size - nProgramOutputs + i], 
                               numbers, ins) + "\n";
    return f;

}
  
//ostream& Print(ostream&) const;
//friend ostream& operator<<(ostream& s, const RbtGPParser &p);

