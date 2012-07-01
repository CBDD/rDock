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

