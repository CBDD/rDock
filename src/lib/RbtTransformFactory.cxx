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


void SetParameterIfExistsInSection(RbtBaseTransform* transform, RbtParameterFileSourcePtr paramsPtr, const RbtString& paramName);
RbtSimAnnTransform* MakeSimmulatedAnnealingTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
RbtGATransform* MakeGeneticAlgorithmTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
RbtAlignTransform* MakeLigandAlignTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
RbtNullTransform* MakeNullTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
RbtRandLigTransform* MakeRandomizeLigandTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
RbtRandPopTransform* MakeRandomizePopulationTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
RbtSimplexTransform* MakeSimplexTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
RbtTransformAgg* MakeAggregateTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
void RegisterScoreFunctionOverridesInTransform(RbtBaseTransform* transform, RbtParameterFileSourcePtr paramsPtr);


// Parameter name which identifies a scoring function definition
RbtString RbtTransformFactory::_TRANSFORM("TRANSFORM");

////////////////////////////////////////
// Constructors/destructors
RbtTransformFactory::RbtTransformFactory() {}

RbtTransformFactory::~RbtTransformFactory() {}

////////////////////////////////////////
// Public methods
////////////////


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

    for (auto& sectionName : transformList) {
        spPrmSource->SetSection(sectionName);
        if (spPrmSource->isParameterPresent(_TRANSFORM)) {
            RbtBaseTransform* transform = MakeTransformFromFile(spPrmSource, sectionName);
            // All examples + the docs instruct to use NullTransform to set the SF overrides, nevertheless
            // it is possible to set the overrides in an arbitrary transform.
            RegisterScoreFunctionOverridesInTransform(transform, spPrmSource);
            pTransformAgg->Add(transform);
        } else if (bThrowError) {
            throw RbtFileMissingParameter(_WHERE_, "Missing " + _TRANSFORM + " parameter in section " + (sectionName));
        }
    }
    return pTransformAgg;
}

// Assumes that the fileSource has the appropriate Section set.
RbtBaseTransform* RbtTransformFactory::MakeTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    RbtString kind = paramsPtr->GetParameterValueAsString(_TRANSFORM);   // Force a cast from RbtVariant to String
    if (kind == RbtSimAnnTransform::_CT) return MakeSimmulatedAnnealingTransformFromFile(paramsPtr, name);
    else if (kind == RbtGATransform::_CT) return MakeGeneticAlgorithmTransformFromFile(paramsPtr, name);
    else if (kind == RbtAlignTransform::_CT) return MakeLigandAlignTransformFromFile(paramsPtr, name);
    else if (kind == RbtNullTransform::_CT) return MakeNullTransformFromFile(paramsPtr, name);
    else if (kind == RbtRandLigTransform::_CT) return MakeRandomizeLigandTransformFromFile(paramsPtr, name);
    else if (kind == RbtRandPopTransform::_CT) return MakeRandomizePopulationTransformFromFile(paramsPtr, name);
    else if (kind == RbtSimplexTransform::_CT) return MakeSimplexTransformFromFile(paramsPtr, name);
    else if (kind == RbtTransformAgg::_CT) return MakeAggregateTransformFromFile(paramsPtr, name);
    else throw RbtBadArgument(_WHERE_, "Unknown transform: " + kind);
}

RbtSimAnnTransform* MakeSimmulatedAnnealingTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    auto transform = new RbtSimAnnTransform(name);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_START_T);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_FINAL_T);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_BLOCK_LENGTH);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_SCALE_CHROM_LENGTH);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_NUM_BLOCKS);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_STEP_SIZE);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_MIN_ACC_RATE);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_PARTITION_DIST);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_PARTITION_FREQ);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimAnnTransform::_HISTORY_FREQ);
    return transform;
}

RbtGATransform* MakeGeneticAlgorithmTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    auto transform = new RbtGATransform(name);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtGATransform::_NEW_FRACTION);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtGATransform::_PCROSSOVER);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtGATransform::_XOVERMUT);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtGATransform::_CMUTATE);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtGATransform::_STEP_SIZE);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtGATransform::_EQUALITY_THRESHOLD);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtGATransform::_NCYCLES);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtGATransform::_NCONVERGENCE);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtGATransform::_HISTORY_FREQ);
    return transform;
}

RbtAlignTransform* MakeLigandAlignTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    auto transform = new RbtAlignTransform(name);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtAlignTransform::_COM);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtAlignTransform::_AXES);
    return transform;
}

// Doesn't have any parameters but let's create it for simmetry;
RbtNullTransform* MakeNullTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    return new RbtNullTransform(name);
}

RbtRandLigTransform* MakeRandomizeLigandTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    auto transform = new RbtRandLigTransform(name);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtRandLigTransform::_TORS_STEP);
    return transform;
}

RbtRandPopTransform* MakeRandomizePopulationTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    auto transform = new RbtRandPopTransform(name);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtRandPopTransform::_POP_SIZE);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtRandPopTransform::_SCALE_CHROM_LENGTH);
    return transform;
}

RbtSimplexTransform* MakeSimplexTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    auto transform = new RbtSimplexTransform(name);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_MAX_CALLS);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_NCYCLES);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_STOPPING_STEP_LENGTH);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_PARTITION_DIST);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_STEP_SIZE);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_CONVERGENCE);
    return transform;
}

RbtTransformAgg* MakeAggregateTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    return new RbtTransformAgg(name);
}

void SetParameterIfExistsInSection(RbtBaseTransform* transform, RbtParameterFileSourcePtr paramsPtr, const RbtString& paramName) {
    if (paramsPtr->isParameterPresent(paramName))
        transform->SetParameter(paramName, paramsPtr->GetParameterValueAsString(paramName));
}

void RegisterScoreFunctionOverridesInTransform(RbtBaseTransform* transform, RbtParameterFileSourcePtr paramsPtr) {
    // Set all the transform parameters from the rest of the parameters listed
    for (auto& paramName : paramsPtr->GetParameterList()) {
        // Look for scoring function request (PARAM@SF). Only SetParamRequest currently supported
        // Parameters of the individual transformers are explicitly set by their respective constructor functions
        // So we only look for score function overrides.
        RbtStringList compList = Rbt::ConvertDelimitedStringToList(paramName, "@");
        if (compList.size() == 2) {
            RbtRequestPtr spReq(new RbtSFSetParamRequest(
                compList[1], compList[0], paramsPtr->GetParameterValueAsString(paramName)
            ));
            transform->AddSFRequest(spReq);
        }
    }
}