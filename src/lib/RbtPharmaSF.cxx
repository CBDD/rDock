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

#include "RbtPharmaSF.h"
#include "RbtWorkSpace.h"
#include "RbtMdlFileSink.h"
#include "RbtLigandError.h"
#include "RbtFileError.h"

//Static data members
RbtString RbtPharmaSF::_CT("RbtPharmaSF");
RbtString RbtPharmaSF::_CONSTRAINTS_FILE("CONSTRAINTS_FILE");
RbtString RbtPharmaSF::_OPTIONAL_FILE("OPTIONAL_FILE");
RbtString RbtPharmaSF::_NOPT("NOPT");
RbtString RbtPharmaSF::_WRITE_ERRORS("WRITE_ERRORS");
    
//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtPharmaSF::RbtPharmaSF(const RbtString& strName) 
  : RbtBaseSF(_CT,strName),m_nopt(0),m_bWriteErrors(false)
{
    //Add parameters It gets the right name in SetupReceptor
    AddParameter(_CONSTRAINTS_FILE,".const");
    AddParameter(_OPTIONAL_FILE,"_opt.const");
    AddParameter(_NOPT,m_nopt);
    AddParameter(_WRITE_ERRORS,m_bWriteErrors);
    SetTrace(1);//Provide a bit of debug output by default
#ifdef _DEBUG
    cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtPharmaSF::~RbtPharmaSF() 
{
#ifdef _DEBUG
    cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
    
void RbtPharmaSF::SetupReceptor() 
{
  m_constrList.clear();
  m_optList.clear();
  m_conScores.clear();
  m_optScores.clear();

  if (GetReceptor().Null())
    return;
  RbtString strWSName = GetWorkSpace()->GetName();
  RbtString strConstraintFile = Rbt::GetRbtFileName("", GetParameter(_CONSTRAINTS_FILE));
  RbtString strOptFile = Rbt::GetRbtFileName("", GetParameter(_OPTIONAL_FILE));
  
  //Create an output sink for ligands that do not have sufficient
  //pharmacophore features. These ligands will be written to an error SD file
  //but will not be docked.
  if (m_bWriteErrors) {
    RbtMolecularFileSinkPtr p = GetWorkSpace()->GetSink();
    if (!p) {
      //This error will be thrown if the receptor model is connected
      //to the workspace BEFORE the ligand output sink
      //Should never happen with the current rbdock.cxx executable
      throw RbtError(_WHERE_, "Ligand output file sink undefined");
    }
    RbtString strOutputFile(p->GetFileName());
    strOutputFile.erase(strOutputFile.find(".sd"),3);
    strOutputFile += "_errors.sd";
    if (GetTrace() > 0) {
      cout << _CT << ": Reading mandatory ph4 constraints from " << strConstraintFile << endl;
      cout << _CT << ": Ligands missing the mandatory features will be written to " << strOutputFile << endl;
      }
    m_spErrorFile = RbtMolecularFileSinkPtr (new RbtMdlFileSink(strOutputFile.c_str(), GetLigand()));
  }
  if (GetTrace() > 0) {
    cout << _CT << ": Reading mandatory ph4 constraints from " << strConstraintFile << endl;
  }
  ifstream constrFile(strConstraintFile.c_str(), ios_base::in);
  if (!constrFile) {
    throw RbtFileReadError(_WHERE_, "cannot open mandatory constraints file " + strConstraintFile);
  }
  Rbt::ZeroCounters();
  Rbt::ReadConstraints(constrFile, m_constrList, true);
  constrFile.close();
  
  //Optional constraints
  ifstream optFile(strOptFile.c_str(), ios_base::in);
  if (optFile.good()) {
    if (GetTrace() > 0) {
      cout << _CT << ": Reading optional ph4 constraints from " << strOptFile << endl;
    }
    Rbt::ReadConstraints(optFile, m_optList, false);
    //Keep m_nopt within range
    SetParameter(_NOPT,std::min(m_nopt,RbtInt(m_optList.size())));
    SetParameter(_NOPT,std::max(m_nopt,0));
    if (GetTrace() > 0) {
      cout << _CT << ": " << m_nopt << " of the optional constraints are required to be satisfied" << endl;
    }
    optFile.close();
  }
  else if (GetTrace() > 0) {
    cout << _CT << ": No optional ph4 constraints file found" << endl;
  }
  
  //Initialise the component score vectors
  m_conScores = RbtDoubleList(m_constrList.size(),0.0);
  m_optScores = RbtDoubleList(m_optList.size(),0.0);
}

void RbtPharmaSF::SetupLigand() {
  if (GetLigand().Null())
    return;
  try {
    if (GetTrace() > 0) {
      cout << _CT << ": Checking mandatory ph4 features..." << endl;
    }
    for (RbtConstraintListIter iter = m_constrList.begin(); iter != m_constrList.end(); iter++) {
      (*iter)->AddAtomList(GetLigand(),true);
    }
    if (GetTrace() > 0) {
      cout << _CT << ": All mandatory features found" << endl;
    }      
  }
  catch (RbtLigandError& e) {
    if (m_bWriteErrors) {
      m_spErrorFile->SetModel(GetLigand());
      m_spErrorFile->Render();
    }
    throw;
    }
  
  //Update the optional features
  for (RbtConstraintListIter iter = m_optList.begin(); iter != m_optList.end(); iter++) {
    (*iter)->AddAtomList(GetLigand(),false);
  }
}

void RbtPharmaSF::SetupScore() 
{
    //No further setup required
}

RbtDouble RbtPharmaSF::RawScore() const 
{
  RbtDouble total(0.0);
  //Store and sum all the mandatory terms
  RbtInt i = 0;
  for (RbtConstraintListConstIter iter = m_constrList.begin(); iter != m_constrList.end(); iter++,i++) {
    m_conScores[i] = (*iter)->Score();
  }
  total = std::accumulate(m_conScores.begin(),m_conScores.end(),total);

  //Store the optional scores for sorting
  i = 0;
  for (RbtConstraintListConstIter iter = m_optList.begin(); iter != m_optList.end(); iter++,i++) {
    m_optScores[i] = (*iter)->Score();
  }
  //partial_sort_copy copies the N lowest scores
  RbtDoubleList lowest(m_nopt);
  std::partial_sort_copy(m_optScores.begin(),m_optScores.end(),lowest.begin(),lowest.end());
  //cout << m_nopt << " lowest optional scores:\t";
  for (RbtInt i = 0; i < lowest.size(); i++) {
    //cout << lowest[i] << "\t";
  }
  //cout << endl;
  total = std::accumulate(lowest.begin(),lowest.end(),total);
  return total;
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtPharmaSF::ParameterUpdated(const RbtString& strName) 
{
  if (strName == _WRITE_ERRORS) {
    m_bWriteErrors = GetParameter(_WRITE_ERRORS);
  }
  else if (strName == _NOPT) {
    m_nopt = GetParameter(_NOPT);
  }
  RbtBaseSF::ParameterUpdated(strName);
}

//Override RbtBaseSF::ScoreMap to provide additional raw descriptors
void RbtPharmaSF::ScoreMap(RbtStringVariantMap& scoreMap) const {
  if (isEnabled()) {
    //Copied from RbtBaseSF
    RbtDouble rs = RawScore();
    RbtString name = GetFullName();
    scoreMap[name] = rs;
    AddToParentMapEntry(scoreMap, rs);
    //Store the mandatory constraint scores
    for (RbtInt i = 0; i < m_conScores.size(); i++) {
      ostrstream field;
      field << name << ".con_" << i+1 << ends;
      scoreMap[field.str()] = m_conScores[i];
      delete field.str();      
    }
    //Store the optional constraint scores (unsorted)
    for (RbtInt i = 0; i < m_optScores.size(); i++) {
      ostrstream field;
      field << name << ".opt_" << i+1 << ends;
      scoreMap[field.str()] = m_optScores[i];
      delete field.str();      
    }
  }
}
