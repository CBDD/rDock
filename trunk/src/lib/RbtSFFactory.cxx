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

#include "RbtSFFactory.h"
//Precalculated-grid scoring functions
#include "RbtVdwGridSF.h"
#include "RbtCavityGridSF.h"
#include "RbtCavityFillSF.h"

//Indexed-grid scoring functions
#include "RbtAromIdxSF.h"
#include "RbtPolarIdxSF.h"
#include "RbtVdwIdxSF.h"

//Intramolecular scoring functions
#include "RbtVdwIntraSF.h"
#include "RbtPolarIntraSF.h"
#include "RbtDihedralIntraSF.h"
#include "RbtDihedralTargetSF.h"


#include "RbtNmrSF.h"

// PMF scoring
#include "RbtSetupPMFSF.h"
#include "RbtPMFIdxSF.h"
#include "RbtPMFGridSF.h"

// SA scoring
#include "RbtSetupSASF.h"
#include "RbtSAIdxSF.h"

// Tethered Atoms scoring
#include "RbtTetherSF.h"

// Pharmacophore scoring
#include "RbtPharmaSF.h"

//Misc scoring functions
#include "RbtConstSF.h"
#include "RbtRotSF.h"

//Pseudo-scoring functions
#include "RbtSetupPolarSF.h"

#include "RbtFileError.h"


//Parameter name which identifies a scoring function definition 
RbtString RbtSFFactory::_SF("SCORING_FUNCTION");

////////////////////////////////////////
//Constructors/destructors
RbtSFFactory::RbtSFFactory() {

}

RbtSFFactory::~RbtSFFactory() {

}

////////////////////////////////////////
//Public methods
////////////////
RbtBaseSF* RbtSFFactory::Create(const RbtString& strSFClass, const RbtString& strName) throw (RbtError) {
  //Precalculated-grid scoring functions
	if (strSFClass == RbtVdwGridSF::_CT) return new RbtVdwGridSF(strName);
	if (strSFClass == RbtCavityGridSF::_CT) return new RbtCavityGridSF(strName);
	if (strSFClass == RbtCavityFillSF::_CT) return new RbtCavityFillSF(strName);

  //Indexed-grid scoring functions
	if (strSFClass == RbtAromIdxSF::_CT) return new RbtAromIdxSF(strName);
	if (strSFClass == RbtPolarIdxSF::_CT) return new RbtPolarIdxSF(strName);
	if (strSFClass == RbtVdwIdxSF::_CT) return new RbtVdwIdxSF(strName);

  //Intramolecular scoring functions
	if (strSFClass == RbtVdwIntraSF::_CT) return new RbtVdwIntraSF(strName);
	if (strSFClass == RbtPolarIntraSF::_CT) return new RbtPolarIntraSF(strName);
	if (strSFClass == RbtDihedralIntraSF::_CT) return new RbtDihedralIntraSF(strName);
	if (strSFClass == RbtDihedralTargetSF::_CT) return new RbtDihedralTargetSF(strName);

	//Restraint scoring functions
	if (strSFClass == RbtNmrSF::_CT) return new RbtNmrSF(strName);
	if (strSFClass == "RbtNmrIntnSF") return new RbtNmrSF(strName);

	// Tethered Atoms scoring
	if (strSFClass == RbtTetherSF::_CT) return new RbtTetherSF(strName);

	// Pharmacophore scoring
	if (strSFClass == RbtPharmaSF::_CT) return new RbtPharmaSF(strName);
	
	// Potential Mean Force (PMF) scoring functions
	if(strSFClass == RbtSetupPMFSF::_CT) return new RbtSetupPMFSF(strName);
	if(strSFClass == RbtPMFIdxSF::_CT) return new RbtPMFIdxSF(strName);
	if(strSFClass == RbtPMFGridSF::_CT) return new RbtPMFGridSF(strName);

	// Solvent Accessible SF
	if(strSFClass == RbtSetupSASF::_CT) return new RbtSetupSASF(strName);
	if(strSFClass == RbtSAIdxSF::_CT) return new RbtSAIdxSF(strName);
	
	//Misc scoring functions
	if (strSFClass == RbtConstSF::_CT) return new RbtConstSF(strName);
	if (strSFClass == RbtRotSF::_CT) return new RbtRotSF(strName);
	
  //Pseudo-scoring functions
	if (strSFClass == RbtSetupPolarSF::_CT) return new RbtSetupPolarSF(strName);

	//Aggregate scoring functions
	if (strSFClass == RbtSFAgg::_CT) return new RbtSFAgg(strName);
	
	throw RbtBadArgument(_WHERE_,"Unknown scoring function "+strSFClass);

}

//Creates an aggregate scoring function from a parameter file source
//Each component SF is in a named section, which should minimally contain a SCORING_FUNCTION parameter
//whose value is the class name to instantiate
//strSFClasses contains a comma-delimited list of SF class names to instantiate
//If strSFClasses is empty, all named sections in spPrmSource are scanned for valid scoring function definitions
//SF parameters are set from the list of parameters in each named section
RbtSFAgg* RbtSFFactory::CreateAggFromFile(RbtParameterFileSourcePtr spPrmSource,
																						const RbtString& strName,
																						const RbtString& strSFClasses) throw (RbtError) {
	//Get list of scoring function objects to create
	RbtStringList sfList = Rbt::ConvertDelimitedStringToList(strSFClasses);
	//If strSFClasses is empty, then default to reading all sections of the
	//parameter file for valid scoring function definitions
	//In this case we do not throw an error if a particular section
	//is not a scoring function, we simply skip it
	RbtBool bThrowError(true);
	if (sfList.empty()) {
		sfList = spPrmSource->GetSectionList();
		bThrowError = false;
	}
	
	//Create empty aggregate
	RbtSFAgg* pSFAgg(new RbtSFAgg(strName));
	
	for (RbtStringListConstIter sfIter = sfList.begin(); sfIter != sfList.end(); sfIter++) {
		spPrmSource->SetSection(*sfIter);
		//Check if this section is a valid scoring function definition
		if (spPrmSource->isParameterPresent(_SF)) {
			RbtString strSFClass(spPrmSource->GetParameterValueAsString(_SF));
			//Create new scoring function according to the string value of _SF parameter
			RbtBaseSF* pSF = Create(strSFClass,*sfIter);
			//Set all the scoring function parameters from the rest of the parameters listed
			RbtStringList prmList = spPrmSource->GetParameterList();
			for (RbtStringListConstIter prmIter = prmList.begin(); prmIter != prmList.end(); prmIter++) {
				if ((*prmIter) != _SF) {//Skip _SF parameter
					pSF->SetParameter(*prmIter,spPrmSource->GetParameterValueAsString(*prmIter));
				}
			}
			pSFAgg->Add(pSF);
		}
		else if (bThrowError) {
			throw RbtFileMissingParameter(_WHERE_,"Missing "+ _SF +" parameter in section "+(*sfIter));
		}
	}
	return pSFAgg;
}

