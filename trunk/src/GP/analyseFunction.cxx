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
