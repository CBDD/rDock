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
// Component transforms
#include "RbtAlignTransform.h"
#include "RbtFileError.h"
#include "RbtGATransform.h"
#include "RbtNullTransform.h"
#include "RbtRandLigTransform.h"
#include "RbtRandPopTransform.h"
#include "RbtSFRequest.h"
#include "RbtSimAnnTransform.h"
#include "RbtSimplexTransform.h"

// Parameter name which identifies a scoring function definition
RbtString RbtTransformFactory::_TRANSFORM("TRANSFORM");

////////////////////////////////////////
// Constructors/destructors
RbtTransformFactory::RbtTransformFactory() {}

RbtTransformFactory::~RbtTransformFactory() {}

////////////////////////////////////////
// Public methods
////////////////

void RbtTransformFactory::AddTransformToAggFromFile(RbtTransformAgg* aggPtr, RbtParameterFileSourcePtr paramsPtr, const RbtString& kind, const RbtString& name) {
    // Create new transform according to the string value of _TRANSFORM parameter
    RbtBaseTransform* pTransform;
    // Component transforms
    if (kind == RbtSimAnnTransform::_CT) pTransform = new RbtSimAnnTransform(name);
    else if (kind == RbtGATransform::_CT) pTransform = new RbtGATransform(name);
    else if (kind == RbtAlignTransform::_CT) pTransform = new RbtAlignTransform(name);
    else if (kind == RbtNullTransform::_CT) pTransform = new RbtNullTransform(name);
    else if (kind == RbtRandLigTransform::_CT) pTransform = new RbtRandLigTransform(name);
    else if (kind == RbtRandPopTransform::_CT) pTransform = new RbtRandPopTransform(name);
    else if (kind == RbtSimplexTransform::_CT) pTransform = new RbtSimplexTransform(name);
    else if (kind == RbtTransformAgg::_CT) pTransform = new RbtTransformAgg(name);
    else throw RbtBadArgument(_WHERE_, "Unknown transform: " + kind);

    // Set all the transform parameters from the rest of the parameters listed
    RbtStringList prmList = paramsPtr->GetParameterList();
    for (RbtStringListConstIter prmIter = prmList.begin(); prmIter != prmList.end(); prmIter++) {
        // Look for scoring function request (PARAM@SF)
        // Only SetParamRequest currently supported
        RbtStringList compList = Rbt::ConvertDelimitedStringToList(*prmIter, "@");
        if (compList.size() == 2) {
            RbtRequestPtr spReq(new RbtSFSetParamRequest(
                compList[1], compList[0], paramsPtr->GetParameterValueAsString(*prmIter)
            ));
            pTransform->AddSFRequest(spReq);
        } else if ((*prmIter) != _TRANSFORM) {  // Skip _TRANSFORM parameter
            pTransform->SetParameter(*prmIter, paramsPtr->GetParameterValueAsString(*prmIter));
        }
    }
    aggPtr->Add(pTransform);
}

// Creates an aggregate transform from a parameter file source
// Each component transform is in a named section, which should minimally contain a TRANSFORM parameter
// whose value is the class name to instantiate
// strTransformClasses contains a comma-delimited list of transform class names to instantiate
// If strTransformClasses is empty, all named sections in spPrmSource are scanned for valid transform definitions
// Transform parameters and scoring function requests are set from the list of parameters in each named section
RbtTransformAgg* RbtTransformFactory::CreateAggFromFile(
    RbtParameterFileSourcePtr spPrmSource, const RbtString& strName, const RbtString& strTransformClasses
) {
    // Get list of transform objects to create
    RbtStringList transformList = Rbt::ConvertDelimitedStringToList(strTransformClasses);
    // If strTransformClasses is empty, then default to reading all sections of the
    // parameter file for valid transform definitions
    // In this case we do not throw an error if a particular section
    // is not a transform, we simply skip it
    RbtBool bThrowError(true);
    if (transformList.empty()) {
        transformList = spPrmSource->GetSectionList();
        bThrowError = false;
    }

    // Create empty aggregate
    RbtTransformAgg* pTransformAgg(new RbtTransformAgg(strName));

    for (RbtStringListConstIter tIter = transformList.begin(); tIter != transformList.end(); tIter++) {
        spPrmSource->SetSection(*tIter);
        // Check if this section is a valid scoring function definition
        if (spPrmSource->isParameterPresent(_TRANSFORM)) {
            RbtString strTransformClass(spPrmSource->GetParameterValueAsString(_TRANSFORM));
            AddTransformToAggFromFile(pTransformAgg, spPrmSource, strTransformClass, *tIter);
        } else if (bThrowError) {
            throw RbtFileMissingParameter(_WHERE_, "Missing " + _TRANSFORM + " parameter in section " + (*tIter));
        }
    }
    return pTransformAgg;
}
