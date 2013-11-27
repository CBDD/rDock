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

#include "RbtTransformFactory.h"
//Component transforms
#include "RbtSimAnnTransform.h"
#include "RbtGATransform.h"
#include "RbtAlignTransform.h"
#include "RbtNullTransform.h"
#include "RbtRandLigTransform.h"
#include "RbtRandPopTransform.h"
#include "RbtSimplexTransform.h"


#include "RbtSFRequest.h"
#include "RbtFileError.h"


//Parameter name which identifies a scoring function definition 
RbtString RbtTransformFactory::_TRANSFORM("TRANSFORM");

////////////////////////////////////////
//Constructors/destructors
RbtTransformFactory::RbtTransformFactory() {

}

RbtTransformFactory::~RbtTransformFactory() {

}

////////////////////////////////////////
//Public methods
////////////////
//Creates a single transform object of type strTransformClass, and name strName
//e.g. strTransformClass = RbtSimAnnTransform
RbtBaseTransform* RbtTransformFactory::Create(const RbtString& strTransformClass, const RbtString& strName) throw (RbtError) {
	//Component transforms
	if (strTransformClass == RbtSimAnnTransform::_CT) return new RbtSimAnnTransform(strName);
	if (strTransformClass == RbtGATransform::_CT) return new RbtGATransform(strName);
	if (strTransformClass == RbtAlignTransform::_CT) return new RbtAlignTransform(strName);
	if (strTransformClass == RbtNullTransform::_CT) return new RbtNullTransform(strName);
	if (strTransformClass == RbtRandLigTransform::_CT) return new RbtRandLigTransform(strName);
	if (strTransformClass == RbtRandPopTransform::_CT) return new RbtRandPopTransform(strName);
	if (strTransformClass == RbtSimplexTransform::_CT) return new RbtSimplexTransform(strName);
	//Aggregate transforms
	if (strTransformClass == RbtTransformAgg::_CT) return new RbtTransformAgg(strName);
	
	throw RbtBadArgument(_WHERE_,"Unknown transform "+strTransformClass);
}

//Creates an aggregate transform from a parameter file source
//Each component transform is in a named section, which should minimally contain a TRANSFORM parameter
//whose value is the class name to instantiate
//strTransformClasses contains a comma-delimited list of transform class names to instantiate
//If strTransformClasses is empty, all named sections in spPrmSource are scanned for valid transform definitions
//Transform parameters and scoring function requests are set from the list of parameters in each named section
RbtTransformAgg* RbtTransformFactory::CreateAggFromFile(RbtParameterFileSourcePtr spPrmSource,
																				const RbtString& strName,
																				const RbtString& strTransformClasses) throw (RbtError) {
	//Get list of transform objects to create
	RbtStringList transformList = Rbt::ConvertDelimitedStringToList(strTransformClasses);
	//If strTransformClasses is empty, then default to reading all sections of the
	//parameter file for valid transform definitions
	//In this case we do not throw an error if a particular section
	//is not a transform, we simply skip it
	RbtBool bThrowError(true);
	if (transformList.empty()) {
		transformList = spPrmSource->GetSectionList();
		bThrowError = false;
	}
	
	//Create empty aggregate
	RbtTransformAgg* pTransformAgg(new RbtTransformAgg(strName));
	
	for (RbtStringListConstIter tIter = transformList.begin(); tIter != transformList.end(); tIter++) {
		spPrmSource->SetSection(*tIter);
		//Check if this section is a valid scoring function definition
		if (spPrmSource->isParameterPresent(_TRANSFORM)) {
			RbtString strTransformClass(spPrmSource->GetParameterValueAsString(_TRANSFORM));
			//Create new transform according to the string value of _TRANSFORM parameter
			RbtBaseTransform* pTransform = Create(strTransformClass,*tIter);
			//Set all the transform parameters from the rest of the parameters listed
			RbtStringList prmList = spPrmSource->GetParameterList();
			for (RbtStringListConstIter prmIter = prmList.begin(); prmIter != prmList.end(); prmIter++) {
				//Look for scoring function request (PARAM@SF)
				//Only SetParamRequest currently supported
				RbtStringList compList = Rbt::ConvertDelimitedStringToList(*prmIter,"@");
				if (compList.size() == 2) {
					RbtRequestPtr spReq(new RbtSFSetParamRequest(compList[1],compList[0],
																											spPrmSource->GetParameterValueAsString(*prmIter)));
					pTransform->AddSFRequest(spReq);
				}
				else if ((*prmIter) != _TRANSFORM) {//Skip _TRANSFORM parameter
					pTransform->SetParameter(*prmIter,spPrmSource->GetParameterValueAsString(*prmIter));
				}
			}
			pTransformAgg->Add(pTransform);
		}
		else if (bThrowError) {
			throw RbtFileMissingParameter(_WHERE_,"Missing "+ _TRANSFORM +" parameter in section "+(*tIter));
		}
	}
	return pTransformAgg;
}
	

