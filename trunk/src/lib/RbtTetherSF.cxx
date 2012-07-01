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

#include "RbtTetherSF.h"
#include "RbtWorkSpace.h"
#include "RbtMdlFileSource.h"
#include <sstream>
using std::istringstream;

//Static data members
RbtString RbtTetherSF::_CT("RbtTetherSF");
RbtString RbtTetherSF::_REFERENCE_FILE("REFERENCE_FILE");
    
//NB - Virtual base class constructor (RbtBaseSF) gets called first,
//implicit constructor for RbtBaseInterSF is called second
RbtTetherSF::RbtTetherSF(const RbtString& strName) 
             : RbtBaseSF(_CT,strName)
{
    //Add parameters It gets the right name in SetupReceptor
    AddParameter(_REFERENCE_FILE,"_reference.sd");
#ifdef _DEBUG
    cout << _CT << " parameterised constructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtTetherSF::~RbtTetherSF() 
{
#ifdef _DEBUG
    cout << _CT << " destructor" << endl;
#endif //_DEBUG
  _RBTOBJECTCOUNTER_DESTR_(_CT);
}
    
RbtIntList RbtTetherSF::ReadTetherAtoms(RbtStringList &strAtoms)
{
    RbtIntList tetherAtoms;
    RbtString strTetherAtoms 
        = Rbt::ConvertListToDelimitedString(strAtoms);
    istringstream ist(strTetherAtoms.c_str());
    RbtInt i;
    char sep[2];
    while (ist >> i)
    {
        tetherAtoms.push_back(i - 1);
        ist.width(2);
        ist >> sep;
    }
    return tetherAtoms;
}

void RbtTetherSF::SetupReceptor() 
{
    if (GetReceptor().Null())
        return;
    RbtString strWSName = GetWorkSpace()->GetName();
    RbtString refExt = GetParameter(_REFERENCE_FILE);
    RbtString refFile = Rbt::GetRbtFileName("", strWSName + refExt);
    RbtMolecularFileSourcePtr spReferenceSD(new RbtMdlFileSource (refFile, false, false, true));
    RbtModelPtr spReferenceMdl(new RbtModel(spReferenceSD));
    RbtStringList strTetherAtomsL = spReferenceMdl->GetDataValue("TETHERED ATOMS");
    RbtIntList tetherAtomsId = ReadTetherAtoms(strTetherAtomsL);
    RbtAtomList refAtoms = spReferenceMdl->GetAtomList();
    for (RbtIntListIter iter = tetherAtomsId.begin(); iter < tetherAtomsId.end() ; iter++) {
      m_tetherCoords.push_back(refAtoms[*iter]->GetCoords());
    }
}

void RbtTetherSF::SetupLigand() 
{
    m_ligAtomList.clear();
    if (GetLigand().Null())
      return;

    RbtStringList strTetherAtomsL = GetLigand()->GetDataValue("TETHERED ATOMS");
    m_tetherAtomList = ReadTetherAtoms(strTetherAtomsL);
    if (m_tetherAtomList.size() != m_tetherCoords.size())
      throw RbtBadArgument(_WHERE_, 
			   "the number of tethered atoms in the ligand SD file should be the same than in the reference SD file");
    
    m_ligAtomList = GetLigand()->GetAtomList();
#ifdef _DEBUG
    cout << _CT << "::SetupLigand(): #ATOMS = " << m_ligAtomList.size() << endl;
#endif //_DEBUG
}

void RbtTetherSF::SetupScore() 
{
  //No further setup required
}

RbtDouble RbtTetherSF::RawScore() const 
{
  RbtDouble score(0.0);
  RbtInt i = 0;
  for (RbtIntListConstIter iter = m_tetherAtomList.begin(); iter < m_tetherAtomList.end() ; iter++, i++)
    score += Rbt::Length2(m_ligAtomList[*iter]->GetCoords(), m_tetherCoords[i]);
  return score;
}

//DM 25 Oct 2000 - track changes to parameter values in local data members
//ParameterUpdated is invoked by RbtParamHandler::SetParameter
void RbtTetherSF::ParameterUpdated(const RbtString& strName) 
{
  RbtBaseSF::ParameterUpdated(strName);
}
