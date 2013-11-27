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
    

