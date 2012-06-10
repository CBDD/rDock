/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/src/GP/analyseFunction.cxx#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../../docs/EULA.txt
***********************************************************************/

#include <stdio.h>
#include "../inc/RbtParser.h"
#include "../inc/RbtGPGenome.h"
#include "../inc/RbtGPPopulation.h"
#include <fstream>
#include <strstream>
void main(int argc,char* argv[])
{
    /*cout << "Genome file: \n";
    RbtString gfile;
    cin >> gfile; */
    RbtString gfile(argv[1]);
    ifstream gstr(gfile.c_str(), ios::in);
    RbtGPGenome gen(gstr);
    /*cout << "Names file: \n";
    RbtString nfile;
    cin >> nfile;*/
    RbtString nfile(argv[2]);
    ifstream nstr(nfile.c_str(), ios::in);
    RbtParser p(gen.GetNIP(), gen.GetNIF(), 
                gen.GetNN(), gen.GetNO());
    RbtInt i = 1;
    RbtInt pos = nstr.tellg();
    while (i!=EOF)
    {
        gen.GetChrom()->Clear();
        nstr.close();
        nstr.open(nfile.c_str(), ios::in);
        cout << p.PrintParse(nstr,gen.GetChrom(),true,false) << endl;
        gen.GetChrom()->Clear();
        nstr.close();
        nstr.open(nfile.c_str(), ios::in);
        cout << p.PrintParse(nstr,gen.GetChrom(),false,false);
        i = EOF;//cin.get();
    }
}
