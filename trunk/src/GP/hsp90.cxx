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

#include <stdio.h>
#include "RbtParser.h"
#include "RbtGPGenome.h"
#include "RbtGPPopulation.h"
#include "RbtGPFitnessFunction.h"
#include "RbtGPFFHSP90.h"
#include "RbtDebug.h"
#include <fstream>
void main(int argc, char* argv[])
{
  try
  {
    RbtRand& theRand = Rbt::GetRbtRand();//ref to random number generator
    if (argc > 1)
        theRand.Seed(atoi(argv[1]));
    else
        theRand.SeedFromClock();
    cout << "Seed: " << theRand.GetSeed() << endl;
    cout << "Output taken from: \n";
    RbtString strInputFile;
    cin >> strInputFile;
    ifstream inputFile(strInputFile.c_str(), ios::in);
    if (!inputFile) 
      throw RbtError(_WHERE_, "can't open " + strInputFile);
    RbtGPGenome::SetStructure(56, 3, 1, 1, 7, 1, 200, 100);
    ifstream desc("descnames", ios::in);
    RbtContextPtr contextp(new RbtCellContext(desc));
    RbtGPFitnessFunctionPtr ff = new RbtGPFFHSP90(contextp);
    cout << "Number of input sets: \n";
    RbtInt nInputSets;
    cin >> nInputSets;
    RbtInt nTestSets = nInputSets / 3;
    RbtString strTrainingFile = "/tmp/beatriz/" + strInputFile + ".training";
    fstream trainingFile(strTrainingFile.c_str(), ios::in|ios::out);
    RbtString strTestFile = "/tmp/beatriz/" + strInputFile + ".test";
    fstream testFile(strTestFile.c_str(), ios::in|ios::out);
    RbtUIntSet randomList;
    cout << nTestSets << endl;
    for (RbtInt i = 0 ; i < nTestSets ; i++)
    {
      RbtBool insertb;
      do
      {
        RbtUInt ran = theRand.GetRandomInt(nInputSets);
        insertb = (randomList.insert(ran)).second;
      }while (!insertb);
    }
    RbtUIntSetIter iter = randomList.end(); 
    for (iter = randomList.begin() ; iter != randomList.end() ; iter++)
      cout << *iter << "\t";
    cout << endl;
    RbtString line;
    getline(inputFile, line);
    testFile << line << endl;
    trainingFile << line << endl;
    for (RbtInt i = 0 ; i < nInputSets ; i++)
    {
      getline(inputFile, line);
      if (randomList.find(i) != randomList.end())
        // found it
        testFile << line << endl;
      else
        trainingFile << line << endl;
    }

    trainingFile.seekp(0, ios::beg);
    testFile.seekp(0, ios::beg);
    
    RbtReturnTypeArray ittrain, sfttrain, ittest, sfttest;
    ff->ReadTables(trainingFile,ittrain,sfttrain);
//    if (!testFile) 
//      throw RbtError(_WHERE_, "can't open" + strTestFile);
    ff->ReadTables(testFile,ittest,sfttest);

    RbtDouble hitlimit = 0.0;
    RbtGPPopulation p(5,5, ff, ittrain, sfttrain);
    p.Initialise(hitlimit,false);
    RbtDouble b = 0.0;
    RbtInt i = 0;
    while ((b < 0.9) && (i < 10000))
    {
        p.EPstep("", 1.0, 0.08, 0.0, 0.0,hitlimit,false);
        b = ff->CalculateFitness(p.Best(),ittrain,sfttrain,false);
        if ((i % 500) == 0)
        {
    //        cout << ff->CalculateFitness(p.Best(),ittest,sfttest,
     //                                    hitlimit,false) << endl;
            cout << b << "\t" 
                 << ff->CalculateFitness(p.Best(),ittrain,sfttrain,
                                         false) << endl;
            if ((i == 5000) && (b < 0))
                break;
        }
        i++;
    }
    cout << i << "\t" << *(p.Best()) << endl;
    cout << "best with training: " << ff->CalculateFitness(p.Best(),ittrain,sfttrain,false) 
                                         << endl;
    cout << "best with test:     " << ff->CalculateFitness(p.Best(),ittest,
                                         sfttest,false) 
                                         << endl;
		inputFile.seekg(0,ios::beg);
    ff->ReadTables(inputFile,ittest,sfttest);
    cout << "best with all:      " << ff->CalculateFitness(p.Best(),ittest,
                                         sfttest,false) 
                                         << endl;

  }
  catch (RbtError& e) 
  {
    cout << e << endl;
  }
  catch (...) 
  {
    cout << "Unknown exception" << endl;
  }

  _RBTOBJECTCOUNTER_DUMP_(cout)
  return 0;

}
