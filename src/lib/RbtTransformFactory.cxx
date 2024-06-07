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


static void SetParameterIfExistsInSection(RbtBaseTransform* transform, RbtParameterFileSourcePtr paramsPtr, const RbtString& paramName);
static RbtSimAnnTransform* MakeSimmulatedAnnealingTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
static RbtGATransform* MakeGeneticAlgorithmTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
static RbtAlignTransform* MakeLigandAlignTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
static RbtNullTransform* MakeNullTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
static RbtRandLigTransform* MakeRandomizeLigandTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
static RbtRandPopTransform* MakeRandomizePopulationTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
static RbtSimplexTransform* MakeSimplexTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
static RbtTransformAgg* MakeAggregateTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name);
static void RegisterScoreFunctionOverridesInTransform(RbtBaseTransform* transform, RbtParameterFileSourcePtr paramsPtr);

static RbtAlignTransform::LigandCenterOfMassPlacementStrategy GetLigandCenterOfMassPlacementStrategyFromFile(RbtParameterFileSourcePtr paramsPtr);
static RbtAlignTransform::LigandAxesAlignmentStrategy GetLigandAxesAlignmentStrategyFromFile(RbtParameterFileSourcePtr paramsPtr);


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
    // If strTransformClasses is empty, then default to reading all sections of the parameter file for valid transform definitions
    // In this case we do not throw an error if a particular section is not a transform, we simply skip it.
    // This is not used anywhere but not sure if any user of the API needs it.
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

static RbtSimAnnTransform* MakeSimmulatedAnnealingTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    const RbtSimAnnTransform::Config& default_config = RbtSimAnnTransform::DEFAULT_CONFIG;
    RbtSimAnnTransform::Config config {
        .initial_temp = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_START_T, default_config.initial_temp),
        .final_temp = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_FINAL_T, default_config.final_temp),
        .num_blocks = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_NUM_BLOCKS, default_config.num_blocks),
        .block_length = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_BLOCK_LENGTH, default_config.block_length),
        .scale_chromosome_length = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_SCALE_CHROM_LENGTH, default_config.scale_chromosome_length),
        .step_size = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_STEP_SIZE, default_config.step_size),
        .min_accuracy_rate = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_MIN_ACC_RATE, default_config.min_accuracy_rate),
        .partition_distance = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_PARTITION_DIST, default_config.partition_distance),
        .partition_frequency = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_PARTITION_FREQ, default_config.partition_frequency),
        .history_frequency = paramsPtr->GetParamOrDefault(RbtSimAnnTransform::_HISTORY_FREQ, default_config.history_frequency),
    };
    return new RbtSimAnnTransform(name, config);
}

static RbtGATransform* MakeGeneticAlgorithmTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    const RbtGATransform::Config& default_config = RbtGATransform::DEFAULT_CONFIG;
    RbtGATransform::Config config {
      .population_size_fraction_as_new_individuals_per_cycle = paramsPtr->GetParamOrDefault(RbtGATransform::_NEW_FRACTION, default_config.population_size_fraction_as_new_individuals_per_cycle),
      .crossover_probability = paramsPtr->GetParamOrDefault(RbtGATransform::_PCROSSOVER, default_config.crossover_probability),
      .cauchy_mutation_after_crossover = paramsPtr->GetParamOrDefault(RbtGATransform::_XOVERMUT, default_config.cauchy_mutation_after_crossover),
      .use_cauchy_distribution_for_mutations = paramsPtr->GetParamOrDefault(RbtGATransform::_CMUTATE, default_config.use_cauchy_distribution_for_mutations),
      .relative_step_size = paramsPtr->GetParamOrDefault(RbtGATransform::_STEP_SIZE, default_config.relative_step_size),
      .equality_threshold = paramsPtr->GetParamOrDefault(RbtGATransform::_EQUALITY_THRESHOLD, default_config.equality_threshold),
      .max_cycles = paramsPtr->GetParamOrDefault(RbtGATransform::_NCYCLES, default_config.max_cycles),
      .num_convergence_cycles = paramsPtr->GetParamOrDefault(RbtGATransform::_NCONVERGENCE, default_config.num_convergence_cycles),
      .history_frequency = paramsPtr->GetParamOrDefault(RbtGATransform::_HISTORY_FREQ, default_config.history_frequency),
    };
    return new RbtGATransform(name, config);
}

static RbtAlignTransform* MakeLigandAlignTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {    
    RbtAlignTransform::Config config {
        .center_of_mass_placement_strategy = GetLigandCenterOfMassPlacementStrategyFromFile(paramsPtr),
        .axes_alignment_strategy = GetLigandAxesAlignmentStrategyFromFile(paramsPtr),
    };
    return new RbtAlignTransform(name, config);
}

// Doesn't have any parameters but let's create it for simmetry;
static RbtNullTransform* MakeNullTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    return new RbtNullTransform(name);
}

static RbtRandLigTransform* MakeRandomizeLigandTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    const RbtRandLigTransform::Config& default_config = RbtRandLigTransform::DEFAULT_CONFIG;
    RbtRandLigTransform::Config config {
        .torsion_step = paramsPtr->GetParamOrDefault(RbtRandLigTransform::_TORS_STEP, default_config.torsion_step),
    };
    return new RbtRandLigTransform(name, config);
}

static RbtRandPopTransform* MakeRandomizePopulationTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    const RbtRandPopTransform::Config& default_config = RbtRandPopTransform::DEFAULT_CONFIG;
    RbtRandPopTransform::Config config {
        .population_size = paramsPtr->GetParamOrDefault(RbtRandPopTransform::_POP_SIZE, default_config.population_size),
        .scale_chromosome_length = paramsPtr->GetParamOrDefault(RbtRandPopTransform::_SCALE_CHROM_LENGTH, default_config.scale_chromosome_length),
    };
    return new RbtRandPopTransform(name, config);
}

static RbtSimplexTransform* MakeSimplexTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    auto transform = new RbtSimplexTransform(name);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_MAX_CALLS);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_NCYCLES);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_STOPPING_STEP_LENGTH);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_PARTITION_DIST);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_STEP_SIZE);
    SetParameterIfExistsInSection(transform, paramsPtr, RbtSimplexTransform::_CONVERGENCE);
    return transform;
}

static RbtTransformAgg* MakeAggregateTransformFromFile(RbtParameterFileSourcePtr paramsPtr, const RbtString& name) {
    return new RbtTransformAgg(name);
}

static void SetParameterIfExistsInSection(RbtBaseTransform* transform, RbtParameterFileSourcePtr paramsPtr, const RbtString& paramName) {
    if (paramsPtr->isParameterPresent(paramName))
        transform->SetParameter(paramName, paramsPtr->GetParameterValueAsString(paramName));
}

static void RegisterScoreFunctionOverridesInTransform(RbtBaseTransform* transform, RbtParameterFileSourcePtr paramsPtr) {
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

static RbtAlignTransform::LigandCenterOfMassPlacementStrategy GetLigandCenterOfMassPlacementStrategyFromFile(RbtParameterFileSourcePtr paramsPtr) {
    if (paramsPtr->isParameterPresent(RbtAlignTransform::_COM)) {
        RbtString placement_strategy_val = paramsPtr->GetParameterValueAsString(RbtAlignTransform::_COM);
        if (placement_strategy_val == "ALIGN")
            return RbtAlignTransform::LigandCenterOfMassPlacementStrategy::COM_ALIGN;
        else if (placement_strategy_val == "RANDOM")
            return RbtAlignTransform::LigandCenterOfMassPlacementStrategy::COM_RANDOM;
        else
            throw RbtBadArgument(_WHERE_, "Invalid ligand center of mass placement strategy: " + placement_strategy_val);
    } else
        return RbtAlignTransform::DEFAULT_CONFIG.center_of_mass_placement_strategy;
}

static RbtAlignTransform::LigandAxesAlignmentStrategy GetLigandAxesAlignmentStrategyFromFile(RbtParameterFileSourcePtr paramsPtr) {
    if (paramsPtr->isParameterPresent(RbtAlignTransform::_AXES)) {
        RbtString alignment_strategy_val = paramsPtr->GetParameterValueAsString(RbtAlignTransform::_AXES);
        if (alignment_strategy_val == "ALIGN")
            return RbtAlignTransform::LigandAxesAlignmentStrategy::AXES_ALIGN;
        else if (alignment_strategy_val == "RANDOM")
            return RbtAlignTransform::LigandAxesAlignmentStrategy::AXES_RANDOM;
        else
            throw RbtBadArgument(_WHERE_, "Invalid ligand axes alignment strategy: " + alignment_strategy_val); 
    } else
        return RbtAlignTransform::DEFAULT_CONFIG.axes_alignment_strategy;
}