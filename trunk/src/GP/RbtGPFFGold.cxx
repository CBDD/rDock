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

#include "RbtGPFFGold.h"
#include "RbtGPGenome.h"
#include "RbtGPParser.h"
#include "RbtDebug.h"
#include <strstream>
#include <fstream>

RbtString RbtGPFFGold::_CT("RbtGPFFGold");

void RbtGPFFGold::ReadTables(istream& in)
{
    RbtReturnType value;
    RbtInt i = 0, j;
    in >> value;
    inputTable.clear();
    SFTable.clear();
    RbtString name;
    while (!in.eof())
    {
        RbtReturnTypeList ivalues;
        ivalues.push_back(new RbtReturnType(value)); 
        for (j = 1 ; j < RbtGPGenome::GetNIP() ; j++)
        {
            in >> value;
            ivalues.push_back(new RbtReturnType(value));
        }
        RbtReturnTypeList sfvalues;
        for (j = 0 ; j < RbtGPGenome::GetNSFI() ; j++)
        {
            in >> value;
            sfvalues.push_back(new RbtReturnType(value));
        }
        //read name, don't store it for now
        in >> name;
//        if (*sfvalues[1] < 3.0)
//        {
            inputTable.push_back(RbtReturnTypeList(ivalues));
            SFTable.push_back(RbtReturnTypeList(sfvalues));
//        }
        i++;
        in >> value;
    }
}

RbtDouble RbtGPFFGold::CalculateFitness(
        RbtGPGenomePtr g, RbtReturnTypeArray& it,
        RbtReturnTypeArray& sft, RbtBool function)
{
    RbtGPParser p(RbtGPGenome::GetNIP(), RbtGPGenome::GetNIF(), 
                RbtGPGenome::GetNN(), RbtGPGenome::GetNO());
    RbtReturnTypeList o;
    RbtDouble tot = 0.0;
    RbtDouble range = 10.0;
    RbtDouble good = 0.0;
    RbtDouble bad = 0.0;
    RbtDouble hitlimit = 2.0;
    for (RbtInt i = 0 ; i < it.size() ; i++)
    {
        RbtReturnTypeList inputs(it[i]);
        RbtReturnTypeList SFValues = sft[i];
        RbtDouble scoreM = *(SFValues[0]);
        RbtDouble rmsd = *(SFValues[1]);
        RbtGPChromosomePtr c(g->GetChrom());
        o = p.Parse(c,inputs);
//        cout << *(o[0]) << endl;
        for (RbtInt j = 0 ; j < RbtGPGenome::GetNO() ; j++)
            if (function)
            {
                RbtDouble d = 2.5 * rmsd + scoreM - *(o[j]);
                tot += d * d ; //abs(d); // * d;
               // if (abs(2.5 * rmsd + scoreM - *(o[j])) < range)
               //     tot++;
            }
            else
                if (*(o[j]) < 0.0) 
                {
                    if (rmsd < hitlimit)
                        good++; // += 1.2;
                    else
                        bad++;
                }
                else
                    if (rmsd >= hitlimit)
                        ;//good++;
                    else
                        good--;
        c->Clear();
    }
        // objective value always an increasing function
    if (function)
        objective = -tot;
    else
//        if (good >=250 && bad <250)
//            objective = 200;
//        else
            objective = good*1-bad;
    // For now, I am using tournament selection. So the fitness
    // function doesn't need to be scaled in any way
    fitness = objective; 
    g->SetFitness(fitness);
    return fitness;
} 
    
RbtDouble RbtGPFFGold::CalculateFitness(
        RbtGPGenomePtr g, RbtReturnTypeArray& it,
        RbtReturnTypeArray& sft, RbtDouble hitlimit, RbtBool function)
{
    RbtGPParser p(g->GetNIP(), g->GetNIF(), g->GetNN(), g->GetNO());
    RbtReturnTypeList o;
    RbtDouble truehits = 0.0;
    RbtDouble falsehits = 0.0;
    RbtDouble truemisses = 0.0;
    RbtDouble falsemisses = 0.0;
    for (RbtInt i = 0 ; i < it.size() ; i++)
    {
        RbtReturnTypeList inputs = it[i];
        RbtReturnTypeList SFValues = sft[i];
        RbtDouble scoreM = *(SFValues[0]);
        RbtDouble rmsd = *(SFValues[1]);
/*        cout << "input  : ";
        for (RbtInt j = 0 ; j < inputs.size() ; j++)
            cout << *(inputs[j]) << " ";
        cout << endl; */
        RbtGPChromosomePtr c = g->GetChrom();
        o = p.Parse(c, inputs);
//        cout << *(o[0]) << endl;
        RbtDouble limit = function ? scoreM + 5 : 0.0;
        for (RbtInt j = 0 ; j < RbtGPGenome::GetNO() ; j++)
        {
            if (rmsd < hitlimit)
                if (*(o[j]) < limit)
                // true hit
                    truehits++;
                else
                    // false miss
                    falsemisses++;
            else
                if (*(o[j]) < limit)
                // false hit
                    falsehits++;
                else
                // true miss
                    truemisses++;
        }
        c->Clear();
    }
        // objective value always an increasing function
    cout << truehits << "\t" << falsehits << "\t" << truemisses 
         << "\t" << falsemisses << endl;
    objective = truehits/(truehits+falsehits);
    fitness = objective;
    g->SetFitness(fitness);
    // For now, I am using tournament selection. So the fitness
    // function doesn't need to be scaled in any way
    return fitness;
}
    

