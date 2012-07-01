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

#include <fstream>

#include "RbtVdwGridSF.h"
#include "RbtFileError.h"

/////////////////////////////////////////////////////////////////////
// MAIN PROGRAM STARTS HERE
/////////////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{
  //Command line arguments and default values
  RbtString strInputFile;
  RbtString strOutputFile("insight.grid");
  RbtInt iGrid=999;

  //Brief help message
  if (argc == 1) {
    cout << endl << "rbconvgrid - converts RbtVdwGridSF binary grid file to InsightII ascii grid file" << endl;
    cout << endl << "Usage:\trbconvgrid -i<InputFile> [-o<OutputFile>] [-n<GridNum>]" << endl;
    cout << endl << "Options:\t-i<InputFile> - input RbtVdwGridSF binary grid filename" << endl;
    cout << "\t\t-o<OutputFile> - output InsightII ascii grid filename" << endl;
    cout << "\t\t-n<GridNum> - grid number to convert (default = list grids)" << endl;
    return 1;
  }

  //Check command line arguments
  cout << endl << "Command line args:" << endl;
  for (RbtInt iarg = 1; iarg < argc; iarg++) {
    cout << argv[iarg];
    RbtString strArg(argv[iarg]);
    if (strArg.find("-i")==0)
      strInputFile = strArg.substr(2);
    else if (strArg.find("-o")==0)
      strOutputFile = strArg.substr(2);
    else if (strArg.find("-n")==0) {
      RbtString strGridNum = strArg.substr(2);
      iGrid = atoi(strGridNum.c_str());
    }
    else {
      cout << " ** INVALID ARGUMENT" << endl;
      return 1;
    }
    cout << endl;
  }

  cout << endl;

  try {
    //Read the grid file header
#if defined(__sgi) && !defined(__GNUC__)
    ifstream istr(strInputFile.c_str(),ios_base::in);
#else
    ifstream istr(strInputFile.c_str(),ios_base::in|ios_base::binary);
#endif
    if (istr) {
      cout << strInputFile << " opened OK" << endl;
    }
    //Read header string
    RbtInt length;
    Rbt::ReadWithThrow(istr, (char*) &length, sizeof(length));
    char* header = new char [length+1];
    Rbt::ReadWithThrow(istr, header, length);
    //Add null character to end of string
    header[length] = '\0';
    //Compare title with 
    RbtBool match = (RbtVdwGridSF::_CT == header);
    delete [] header;
    if (!match) {
      throw RbtFileParseError(_WHERE_,"Invalid title string in " + strInputFile);
    }
    
    //Skip the appropriate number of grids
    RbtInt nGrids;
    Rbt::ReadWithThrow(istr, (char*) &nGrids, sizeof(nGrids));
    cout << "File contains " << nGrids << " grids..." << endl;
    if ((iGrid > nGrids) || (iGrid < 1)) {
      cout << "Listing grids..." << endl;
    }
    else {
      cout << "Locating grid# " << iGrid << "..." << endl;
    }
    RbtRealGridPtr spGrid;
    for (RbtInt i = 1; (i <= nGrids) && (i <= iGrid); i++) {
      //Read the atom type string
      Rbt::ReadWithThrow(istr, (char*) &length, sizeof(length));
      char* szType = new char [length+1];
      Rbt::ReadWithThrow(istr, szType, length);
      //Add null character to end of string
      szType[length] = '\0';
      RbtString strType(szType);
      delete [] szType;
      RbtTriposAtomType triposType;
      RbtTriposAtomType::eType aType = triposType.Str2Type(strType);
      cout << "Grid# " << i << "\t" << "atom type=" << strType << " (type #" << aType << ")" << endl;
      spGrid = RbtRealGridPtr(new RbtRealGrid(istr));
    }
    istr.close();
    //If we are not in listing mode, write the grid
    if ((iGrid <= nGrids) && (iGrid >= 1)) {
      cout << "Writing grid# " << iGrid << " to " << strOutputFile << "..." << endl;
      ofstream ostr(strOutputFile.c_str());
      spGrid->PrintInsightGrid(ostr);
      ostr.close();
    }
  }
  catch (RbtError& e) {
    cout << e << endl;
  }
  catch (...) {
    cout << "Unknown exception" << endl;
  }
  
  _RBTOBJECTCOUNTER_DUMP_(cout)
    
  return 0;
}
