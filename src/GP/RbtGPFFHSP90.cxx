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

#include "RbtGPFFHSP90.h"
#include "RbtGPGenome.h"
#include "RbtGPParser.h"
#include "RbtParser.h"
#include "RbtDebug.h"
#include "RbtFilterExpression.h"
#include "RbtFilterExpressionVisitor.h"
#include "RbtCellTokenIter.h"
#include "RbtTokenIter.h"
#include "RbtRand.h"
#include <strstream>
#include <fstream>
#include <cassert>

RbtString RbtGPFFHSP90::_CT("RbtGPFFHSP90");
int nInversions(RbtInt idx, RbtReturnTypeArray& list);

void RbtGPFFHSP90::ReadTables(istream& in, RbtReturnTypeArray& it,
                              RbtReturnTypeArray& sft)
{
    RbtReturnType value;
    RbtInt nip, nsfi;
    RbtInt i = 0, j;
    in >> nip ;
    in.get();
    RbtInt nctes = 15;
    RbtGPGenome::SetNIP(nip+nctes);
    in >> nsfi;
    RbtGPGenome::SetNSFI(nsfi);
    inputTable.clear();
    SFTable.clear();
    char name[80];
        //read mol, don't store it for now
    in.getline(name, 80, ',');
    CreateRandomCtes(nctes);
    while (!in.eof())
    {
        RbtReturnTypeList sfvalues;
        for (j = 0 ; j < nsfi ; j++)
        {
            in >> value;
            in.get(); // read ','
            sfvalues.push_back(new RbtReturnType(value));
        }
        //read Ki, don't store it
        in.getline(name, 80, ',');
        //read Name, don't store it
        in.getline(name, 80, ',');
        //read RMSD, don't store it
        in.getline(name, 80, ',');
        in >> value;
        RbtReturnTypeList ivalues;
        ivalues.push_back(new RbtReturnType(value)); 
        for (j = 1 ; j < nip ; j++)
        {
            in.get(); // read ','
            in >> value;
            ivalues.push_back(new RbtReturnType(value));
        }
        for (j = 0; j < nctes ; j++)
        {
            ivalues.push_back(new RbtReturnType(ctes[j]));
        }
        inputTable.push_back(RbtReturnTypeList(ivalues));
        SFTable.push_back(RbtReturnTypeList(sfvalues));
        i++;
        in.getline(name, 80, ',');
    }
    it = inputTable;
    sft = SFTable;

}

RbtDouble RbtGPFFHSP90::CalculateFitness(RbtGPGenomePtr g, 
                                         RbtReturnTypeArray& it,
                                         RbtReturnTypeArray& sft,
//                                         RbtDouble limit,
																				 RbtBool function)
{
    RbtReturnTypeArray o;
//    o.push_back(new RbtReturnType(1.1));
    RbtReturnType oldo;
    RbtDouble tot = 0.0;
    RbtDouble good = 0.0;
    RbtDouble bad = 0.0;
    RbtGPChromosomePtr c(g->GetChrom());
    RbtParser p2;
    RbtCellTokenIterPtr ti(new RbtCellTokenIter(c, contextp));
    RbtFilterExpressionPtr fe(p2.Parse(ti, contextp));
    RbtDouble meanReal = 0.0, meanPred = 0.0;
    for (RbtInt i = 0 ; i < it.size() ; i++)
    {
        RbtReturnTypeList inputs(it[i]);
        RbtReturnTypeList pred;
        pred.clear();
        for (RbtInt j = 0 ; j < inputs.size() ; j++)
            contextp->Assign(j, *(inputs[j]));
        RbtReturnTypeList SFValues = sft[i];
        RbtDouble hit = *(SFValues[SFValues.size()-1]);
        EvaluateVisitor visitor(contextp);
        fe->Accept(visitor);
        pred.push_back(new RbtReturnType(fe->GetValue()));
        o.push_back(pred);
//        cout << *(SFValues[0]) << "\t" << *(o[i][0]) << endl;
    }
		RbtInt ninv = nInversions(0, o);
//		cout << "n inversions " << ninv << endl;
/*    RbtDouble sumDiff = 0.0, 
              sumSqrReal = 0.0, sumSqrPred = 0.0;
    
    meanReal = meanReal / it.size();
    meanPred = meanPred / it.size();
//    cout << endl;
//    cout << it.size() << "\t" << meanPred << "\t" << meanReal << endl;
    for (RbtInt i = 0 ; i < it.size() ; i++)
    {
      RbtDouble t1 = *sft[i][0]; // - meanReal;
      RbtDouble t2 = *o[i][0]; // - meanPred;
      sumDiff += abs(t1 - t2);
//    cout << sumDiff << "\t" << sumSqrReal << "\t" << sumSqrPred << endl;
    }
    objective = sumDiff / it.size();
                    
//    cout << objective << endl; */

    // For now, I am using tournament selection. So the fitness
    // function doesn't need to be scaled in any way
        // objective value increasing function
    fitness = -ninv; 
    g->SetFitness(fitness);
    return fitness;
} 
    
RbtDouble RbtGPFFHSP90::CalculateFitness(
        RbtGPGenomePtr g, RbtReturnTypeArray& it,
        RbtReturnTypeArray& sft, RbtDouble limit,RbtBool function)
{
    RbtReturnTypeArray o;
    RbtReturnType oldo;
    RbtDouble tot = 0.0;
    RbtDouble good = 0.0;
    RbtDouble bad = 0.0;
    RbtGPChromosomePtr c(g->GetChrom());
    RbtParser p2;
    RbtCellTokenIterPtr ti(new RbtCellTokenIter(c, contextp));
    RbtFilterExpressionPtr fe(p2.Parse(ti, contextp));
    RbtDouble meanReal = 0.0, meanPred = 0.0;
    for (RbtInt i = 0 ; i < it.size() ; i++)
    {
        RbtReturnTypeList inputs(it[i]);
        RbtReturnTypeList pred;
        pred.clear();
        for (RbtInt j = 0 ; j < inputs.size() ; j++)
            contextp->Assign(j, *(inputs[j]));
        RbtReturnTypeList SFValues = sft[i];
        RbtDouble hit = *(SFValues[SFValues.size()-1]);
        EvaluateVisitor visitor(contextp);
        fe->Accept(visitor);
        pred.push_back(new RbtReturnType(fe->GetValue()));
        o.push_back(pred);
        meanPred += *(o[i][0]);
        meanReal += *(SFValues[0]);
        cout << *(SFValues[0]) << "\t" << *(o[i][0]) << endl;
    }
    RbtDouble sumDiff = 0.0, 
              sumSqrReal = 0.0, sumSqrPred = 0.0;
    
    meanReal = meanReal / it.size();
    meanPred = meanPred / it.size();
//    cout << endl;
    for (RbtInt i = 0 ; i < it.size() ; i++)
    {
      RbtDouble t1 = *sft[i][0] - meanReal;
      RbtDouble t2 = *o[i][0] - meanPred;
      sumDiff += t1 * t2;
      sumSqrReal += t1 * t1;
      sumSqrPred += t2 * t2;
    }
    objective = (sumDiff * sumDiff) 
                    / (sumSqrReal * sumSqrPred) ;
    // For now, I am using tournament selection. So the fitness
    // function doesn't need to be scaled in any way
        // objective value increasing function
    fitness = objective; 
    g->SetFitness(fitness);
    return fitness;
}
    
void RbtGPFFHSP90::CreateRandomCtes(RbtInt nctes)
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

int nInversions(RbtInt idx, RbtReturnTypeArray& list)
{
	RbtInt count = 0;
	RbtInt len = list.size();
	if (idx == (len - 1))
		return 0;
	for (RbtInt i = idx + 1 ; i < len ; i++)
	{
		if (*(list[idx][0]) >= *(list[i][0]))
			count++;
	}
	return (count + nInversions(idx+1, list));
}
