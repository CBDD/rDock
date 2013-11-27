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

#include "RbtGPFFCrossDock.h"
#include "RbtGPGenome.h"
#include "RbtGPParser.h"
#include "RbtDebug.h"
#include <strstream>
#include <fstream>

RbtString RbtGPFFCrossDock::_CT("RbtGPFFCrossDock");

void RbtGPFFCrossDock::ReadTables(istream& in, RbtReturnTypeArray& it,
                                           RbtReturnTypeArray& sft)
{
    RbtReturnType value;
    RbtInt nip, nsfi;
    RbtInt i = 0, j, recordn;
    in >> nip ;
    in.get();
    RbtInt nctes = 15;
    RbtGPGenome::SetNIP(nip+nctes);
    in >> nsfi;
    RbtGPGenome::SetNSFI(nsfi);
    in >> recordn;
    inputTable.clear();
    SFTable.clear();
    char name[80];
    CreateRandomCtes(nctes);
    while (!in.eof())
    {
        //read name, don't store it for now
        in.get();
        in.getline(name, 80, ',');
        in >> value;
        RbtReturnTypeList ivalues;
        ivalues.push_back(new RbtReturnType(value)); 
        for (j = 1 ; j < nip ; j++)
        {
            in.get();
            in >> value;
            ivalues.push_back(new RbtReturnType(value));
        }
        for (j = 0; j < nctes ; j++)
            ivalues.push_back(new RbtReturnType(ctes[j]));
        RbtReturnTypeList sfvalues;
        for (j = 0 ; j < nsfi ; j++)
        {
            in.get();
            in >> value;
            sfvalues.push_back(new RbtReturnType(value));
        }
        inputTable.push_back(RbtReturnTypeList(ivalues));
        SFTable.push_back(RbtReturnTypeList(sfvalues));
        i++;
        in >> recordn;
    }
    cout << "Read: " << inputTable[0][0] << endl;
    it = inputTable;
    sft = SFTable;
}

RbtDouble RbtGPFFCrossDock::CalculateFitness(RbtGPGenomePtr g, 
                                         RbtReturnTypeArray& it,
                                         RbtReturnTypeArray& sft,
                                         RbtBool function)
{
    if (function)
    {
        cout << "Error, no function possible with Cross Docking\n";
        exit(1);
    }
    RbtGPParser p(g->GetNIP(), g->GetNIF(), g->GetNN(), g->GetNO());
    RbtReturnTypeList o;
    RbtDouble good = 0.0;
    RbtDouble bad = 0.0;
    RbtDouble neutral = 0.0;
    RbtDouble hitlimit = 0.0;
    RbtInt ntm, nm;
    RbtReturnTypeList inputs;
    RbtReturnTypeList SFValues;
    RbtDouble hit;
    RbtInt i = 0;
    while (i < it.size())
    {
        inputs = it[i];
        SFValues = sft[i];
        hit = *(SFValues[SFValues.size()-1]);
        ntm = 0;
        nm = 0;
        while (hit > hitlimit)
        {
/*            cout << "input  : ";
            for (RbtInt j = 0 ; j < inputs.size() ; j++)
                cout << *(inputs[j]) << " ";
            cout << endl; */
            RbtGPChromosomePtr c = g->GetChrom();
            o = p.Parse(c, inputs);
            c->Clear();
            //cout << *(o[0]) << endl;
            if (*(o[0]) >= hitlimit)
                ntm++;
            nm++;
            i++;
            inputs = it[i];
            SFValues = sft[i];
            hit = *(SFValues[SFValues.size()-1]);
        }
        RbtGPChromosomePtr c = g->GetChrom();
        o = p.Parse(c, inputs);
        c->Clear();
        if (*(o[0]) < hitlimit)
        {
     //       if (ntm >= (nm / 2.0))
            if (ntm >= ((nm * 2.0) / 3.0))
                good++;
            else
                neutral++;
        }
        else
		{
    //        if (ntm >= (nm / 2.0))
    //        if (ntm >= ((nm * 2.0) / 3.0))
                bad++;
		}
        i++;
    }
    objective = good-0.4*neutral-bad;
    fitness = objective;
    g->SetFitness(fitness);
    // For now, I am using tournament selection. So the fitness
    // function doesn't need to be scaled in any way
    return fitness;
}
    
RbtDouble RbtGPFFCrossDock::CalculateFitness(
        RbtGPGenomePtr g, RbtReturnTypeArray& it,
        RbtReturnTypeArray& sft, RbtDouble hitlimit,RbtBool function)
{
    RbtGPParser p(g->GetNIP(), g->GetNIF(), g->GetNN(), g->GetNO());
    RbtReturnTypeList o;
    RbtDouble good = 0.0;
    RbtDouble bad = 0.0;
    RbtDouble neutral = 0.0;
    RbtInt ntm, nm;
    RbtReturnTypeList inputs;
    RbtReturnTypeList SFValues;
    RbtDouble hit;
    RbtInt i = 0;
    while (i < it.size())
    {
        inputs = it[i];
        SFValues = sft[i];
        hit = *(SFValues[SFValues.size()-1]);
        ntm = 0;
        nm = 0;
        while (hit > hitlimit)
        {
/*            cout << "input  : ";
            for (RbtInt j = 0 ; j < inputs.size() ; j++)
                cout << *(inputs[j]) << " ";
            cout << endl; */
            RbtGPChromosomePtr c = g->GetChrom();
            o = p.Parse(c, inputs);
            c->Clear();
            //cout << *(o[0]) << endl;
            RbtDouble limit = 0.0;
            if (*(o[0]) >= limit)
                ntm++;
            nm++;
            i++;
            inputs = it[i];
            SFValues = sft[i];
            hit = *(SFValues[SFValues.size()-1]);
        }
        RbtGPChromosomePtr c = g->GetChrom();
        o = p.Parse(c, inputs);
        c->Clear();
        if (*(o[0]) < hitlimit)
        {
            if (ntm >= ((nm * 2.0) / 3.0))
                good++;
            else
                neutral++;
        }
        else
            bad++;
        i++;
    }
        // objective value always an increasing function
    cout << good << "\t" << neutral << "\t" << bad << endl;
    objective = good/(good+bad);
    fitness = objective;
    //g->SetFitness(fitness);
    // For now, I am using tournament selection. So the fitness
    // function doesn't need to be scaled in any way
    return fitness;
}
    
void RbtGPFFCrossDock::CreateRandomCtes(RbtInt nctes)
{
    if (ctes.size() == 0) // it hasn't already been initialized
    {
        RbtInt a, b;
        RbtDouble c;
        ctes.push_back(0.0);
        ctes.push_back(1.0);
        cout << "c0 \t0.0" << endl;
        cout << "c1 \t1.0" << endl;
        for (RbtInt i = 0 ; i < (nctes-2) ; i++)
        {
            a = m_rand.GetRandomInt(200) - 100;
            b = m_rand.GetRandomInt(10) - 5;
            c = (a/10.0)*pow(10,b);
            cout << "c" << i+2 << " \t" << c << endl;
            ctes.push_back(c);
        }
    }
}

