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
