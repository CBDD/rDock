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

//Calculates atomic RMSD of each SD record with reference structure
//Requires Daylight SMARTS toolkit to check all symmetry-related atomic numbering paths
//Reference structure can be a substructural fragment of each ligand in the SD file if desired
//No requirement for consistent numbering scheme between reference structure and SD file

#include <sstream>

#include "RbtMdlFileSource.h"
#include "RbtMdlFileSink.h"
#include "RbtModelError.h"
#include "RbtModel.h"
#include "RbtSmarts.h"

typedef vector<RbtCoordList> RbtCoordListList;
typedef RbtCoordListList::iterator RbtCoordListListIter;
typedef RbtCoordListList::const_iterator RbtCoordListListConstIter;

//RMSD calculation between two coordinate lists
RbtDouble rmsd(const RbtCoordList& rc, const RbtCoordList& c) {
  RbtInt nCoords = rc.size();
  if (c.size() != nCoords) {
    return 0.0;
  }
  else {
    RbtDouble rms(0.0);
    for (RbtInt i = 0; i < nCoords; i++) {
      rms += Rbt::Length2(rc[i],c[i]);
    }
    rms = sqrt(rms/float(nCoords));
    return rms;
  }
}

int main(int argc,char* argv[])
{
  if (argc < 3) {
    cout << "smart_rms <ref sdfile> <input sdfile> [<output sdfile>]" << endl;
    cout << "RMSD is calculated for each record in <input sdfile> against <ref sdfile> (heavy atoms only)" << endl;
    cout << "If <output sdfile> is defined, records are written to output file with RMSD data field" << endl;
    cout << endl << "NOTE:" << endl;
    cout << "\tThis version uses the Daylight SMARTS toolkit to check all symmetry-related atom numbering paths" << endl;
    cout << "\tto ensure that the true RMSD is reported. Will only run on licensed machines (giles)" << endl;
    return 1;
  }

  RbtString strRefSDFile(argv[1]);
  RbtString strInputSDFile(argv[2]);
  RbtString strOutputSDFile;
  RbtBool bOutput(false);
  if (argc > 3) {
    strOutputSDFile = argv[3];
    bOutput = true;
  }
  std::ios_base::fmtflags oldflags = cout.flags();//save state

  RbtDoubleList scoreVec;
  RbtDoubleList rmsVec;

  try {
    RbtMolecularFileSourcePtr spRefFileSource(new RbtMdlFileSource(Rbt::GetRbtFileName("data/ligands",strRefSDFile),false,false,true));
    //DM 16 June 2006 - remove any solvent fragments from reference
    //The largest fragment in each SD record always has segment name="H"
    //for reasons lost in the mists of RiboDock history
    spRefFileSource->SetSegmentFilterMap(Rbt::ConvertStringToSegmentMap("H"));
    //Get reference ligand (first record)
    RbtModelPtr spRefModel(new RbtModel(spRefFileSource));
    RbtCoordListList cll;
    RbtString strSmarts;
    RbtString strSmiles;
    RbtAtomListList pathset = DT::QueryModel(spRefModel, strSmarts, strSmiles);
    cout << "Reference SMILES: " << strSmiles << endl;
    cout << "Paths found = " << pathset.size() << endl;
    if (pathset.empty()) {
      exit(0);
    }
    //Use the SMILES string for the reference to query each record in the SD file
    //This has the useful side effect that the numbering scheme in the SD file
    //does not have to be consistent with that in the reference structure
    //Also allows the reference to be a substructural fragment of each ligand
    strSmarts = strSmiles;
    for (RbtAtomListListConstIter iter = pathset.begin(); iter != pathset.end(); iter++) {
      RbtCoordList coords;
      Rbt::GetCoordList(*iter,coords);
      cll.push_back(coords);
    }
    RbtInt nCoords = cll.front().size();

    cout << "molv_	rms rms	rmc rmc" << endl;//Dummy header line to be like do_anal

    cout.precision(3);
    cout.setf(ios_base::fixed,ios_base::floatfield);

    ///////////////////////////////////
    //MAIN LOOP OVER LIGAND RECORDS
    ///////////////////////////////////
    RbtMolecularFileSourcePtr spMdlFileSource(new RbtMdlFileSource(Rbt::GetRbtFileName("data/ligands",strInputSDFile),false,false,true));
    RbtMolecularFileSinkPtr spMdlFileSink;
    if (bOutput) {
      spMdlFileSink = new RbtMdlFileSink(strOutputSDFile,RbtModelPtr());
    }
    for (RbtInt nRec=1; spMdlFileSource->FileStatusOK(); spMdlFileSource->NextRecord(), nRec++) {
      RbtError molStatus = spMdlFileSource->Status();
      if (!molStatus.isOK()) {
	cout << molStatus << endl;
	continue;
      }
      //DM 16 June 2006 - remove any solvent fragments from each record
      spMdlFileSource->SetSegmentFilterMap(Rbt::ConvertStringToSegmentMap("H"));
      RbtModelPtr spModel(new RbtModel(spMdlFileSource));
      RbtAtomListList pathset1 = DT::QueryModel(spModel, strSmarts, strSmiles);
      //cout << "SMILES: " << strSmiles << endl;
      //cout << "Paths found = " << pathset1.size() << endl;
      if (pathset1.empty()) {
	continue;
      }
      //we only need to retrieve the coordinate list for the first matching atom numbering path
      //as we have already stored all the alternative numbering schemes for the reference structure
      RbtCoordList coords;
      Rbt::GetCoordList(pathset1.front(),coords);
	
      if (coords.size() == nCoords) {//Only calculate RMSD if atom count is same as reference
	RbtDouble rms(9999.9);
	for (RbtCoordListListConstIter cIter = cll.begin(); cIter != cll.end(); cIter++) {
	  RbtDouble rms1 = rmsd(*cIter,coords);
	  //cout << "\tRMSD = " << rms1 << endl;
	  rms = std::min(rms,rms1);
	}
	spModel->SetDataValue("RMSD",rms);
	RbtDouble score = spModel->GetDataValue("SCORE");
	RbtDouble scoreInter = spModel->GetDataValue("SCORE.INTER");
	RbtDouble scoreIntra = spModel->GetDataValue("SCORE.INTRA");
	scoreVec.push_back(score);
	rmsVec.push_back(rms);

	cout << nRec << "\t" << score << "\t" << scoreInter << "\t" << scoreIntra << "\t" << rms << "\t" << 0.0 << endl; 
	if (bOutput) {
	  spMdlFileSink->SetModel(spModel);
	  spMdlFileSink->Render();
	}
      }
    }
    //END OF MAIN LOOP OVER LIGAND RECORDS
    ////////////////////////////////////////////////////
    RbtDoubleListConstIter sIter = scoreVec.begin();
    RbtDoubleListConstIter rIter = rmsVec.begin();
    RbtDouble minScore = *std::min_element(scoreVec.begin(),scoreVec.end());
    RbtDouble zGood(0.0);
    RbtDouble zPartial(0.0);
    RbtDouble zBad(0.0);
    for (; (sIter != scoreVec.end()) && (rIter != rmsVec.end()); sIter++, rIter++) {
      RbtDouble de = (*sIter) - minScore;
      RbtDouble z = exp(-de/(8.314e-3 * 298.0));
      RbtDouble r = (*rIter);
      if (r < 2.05)
	zGood += z;
      else if (r < 3.05)
	zPartial += z;
      else
	zBad += z;
    }
    //cout << "zGood,zPartial,zBad = " << zGood << "," << zPartial << "," << zBad << endl;
  }
  catch (RbtError& e) {
    cout << e << endl;
  }
  catch (...) {
    cout << "Unknown exception" << endl;
  }

  cout.flags(oldflags);//Restore state

  _RBTOBJECTCOUNTER_DUMP_(cout)

  return 0;
}
