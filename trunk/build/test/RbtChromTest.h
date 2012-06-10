//Unit tests for RbtChromElement and subclasses
//
//Required input files:
//1koc.prm  rDock receptor file
//1koc.psf  Receptor topology file
//1koc.crd  Receptor coordinate file
//1koc_c.sd Ligand coordinate file
//1koc.as   Docking site
//
//Required environment:
//Make sure the above files are colocated in a single directory
//and define RBT_HOME env. variable to point at this directory
#ifndef RBTCHROMTEST_H_
#define RBTCHROMTEST_H_

#include <cppunit/extensions/HelperMacros.h>

#include "RbtChrom.h"
#include "RbtModel.h"
#include "RbtDockingSite.h"
#include "RbtBiMolWorkSpace.h"
#include "RbtSFAgg.h"

class RbtChromTest : public CppUnit::TestFixture  {
CPPUNIT_TEST_SUITE( RbtChromTest );
CPPUNIT_TEST( testChromFactory );
CPPUNIT_TEST( testCloneDestructor );
CPPUNIT_TEST( testCloneLengthEquals );
CPPUNIT_TEST( testChromOperatorEquals );
CPPUNIT_TEST( testChromOperatorEqualsAfterMutate );
CPPUNIT_TEST( testSyncToModel );
CPPUNIT_TEST( testSyncToModelAfterMutate );
CPPUNIT_TEST( testSyncToModelAfterReset );
CPPUNIT_TEST( testCloneReset );
CPPUNIT_TEST( testMutatedCloneReset );
CPPUNIT_TEST( testChromOperatorEqualsAfterRandomise );
CPPUNIT_TEST( testSyncToModelAfterRandomise );
CPPUNIT_TEST( testSyncToModelAfterRandomiseReset );
CPPUNIT_TEST( testChromOperatorEqualsAfterMutateReset );
CPPUNIT_TEST( testChromOperatorEqualsAfterRandomiseReset );
CPPUNIT_TEST( testNullCrossover );
CPPUNIT_TEST( testEuler );
CPPUNIT_TEST( testRepeatedSync );
CPPUNIT_TEST( testCompareWithNullChrom );
CPPUNIT_TEST( testCompareAfterMutate );
CPPUNIT_TEST( testOperatorEqualsWithModifiedThreshold );
CPPUNIT_TEST( testPopulationConstructor );
CPPUNIT_TEST_EXCEPTION( testPopulationConstructorZeroSize, RbtBadArgument );
CPPUNIT_TEST_EXCEPTION( testPopulationConstructorNullChromosome, RbtBadArgument );
CPPUNIT_TEST_EXCEPTION( testPopulationConstructorNullSF, RbtBadArgument );
CPPUNIT_TEST( testPopulationRWFitness );
CPPUNIT_TEST( testPopulationRWSelect );
CPPUNIT_TEST( testPopulationGAstep );
CPPUNIT_TEST( testRbtModelGetChromUndefinedFlexData );
CPPUNIT_TEST( testRbtModelGetChromNullFlexData );
CPPUNIT_TEST( testRbtModelGetChrom );
CPPUNIT_TEST( testRandomiseTetheredOrientation );
CPPUNIT_TEST( testRandomiseTetheredCOM );
CPPUNIT_TEST( testMutateTetheredOrientation );
CPPUNIT_TEST( testMutateTetheredCOM );
CPPUNIT_TEST( testCrossoverTetheredOrientation );
CPPUNIT_TEST( testCrossoverTetheredCOM );
CPPUNIT_TEST( testRandomiseTetheredDihedral );
CPPUNIT_TEST( testMutateTetheredDihedral );
CPPUNIT_TEST( testCrossoverTetheredDihedral );
CPPUNIT_TEST( testRandomiseOccupancy );
CPPUNIT_TEST( testOccupancyThreshold );
CPPUNIT_TEST_SUITE_END();

public:
  static const RbtDouble TINY = 1E-4;
  //TextFixture methods
  void setUp();
  void tearDown();
  
  //extra setup for building a simple scoring function and workspace
  void setupWorkSpace();

  //rdock helper methods
  //RMSD calculation between two coordinate lists
  RbtDouble rmsd(const RbtCoordList& rc, const RbtCoordList& c);
  //Measures the variation in chromosome vector values for repeated
  //mutations or randomisations
  void measureRandOrMutateDiff(RbtChromElement* chrom,
                            RbtInt nTrials,
                            RbtBool bMutate,//true=mutate; false=randomise
                            RbtDouble& meanDiff,
                            RbtDouble& minDiff,
                            RbtDouble& maxDiff);
  //Measures the variation in chromosome vector values for repeated
  //crossovers
  void measureCrossoverDiff(RbtChromElement* chrom,
                            RbtInt nTrials,
                            RbtDouble& meanDiff,
                            RbtDouble& minDiff,
                            RbtDouble& maxDiff);
                            
 
  //The unit tests
  //1) Does RbtChromFactory create a valid chromosome?
  void testChromFactory();
  //2) Can we delete a cloned chromosome without affecting the original?
  void testCloneDestructor();
  //3) Is the length of a cloned chromosome equal to the original?
  void testCloneLengthEquals();
  //4) Does operator== return true for a clone and an original (before mutation)?
  void testChromOperatorEquals();
  //5) Does operator!= return true for a mutated clone and an original?
  void testChromOperatorEqualsAfterMutate();
  //6) Does SyncToModel leave the model coords unchanged if chromosome has
  //not been mutated?
  void testSyncToModel();
  //7) Does SyncToModel change the model coords after a mutation?
  void testSyncToModelAfterMutate();
  //8) Does SyncToModel leave the model coords unchanged if chromosome has
  //been mutated then reset?
  void testSyncToModelAfterReset();
  //9) If we make a clone BEFORE the original has been mutated, check that
  //Reset() on the clone will reset to the original (unmutated) model coords 
  void testCloneReset();
  //10) If we make a clone AFTER the original has been mutated, check that
  //Reset() on the clone still resets to the original (unmutated) model coords 
  void testMutatedCloneReset();
  //11) Does operator!= return true for a randomised clone and an original?
  void testChromOperatorEqualsAfterRandomise();
  //12) Does SyncToModel change the model coords after a randomisation?
  void testSyncToModelAfterRandomise();
  //13) Does SyncToModel leave the model coords unchanged if chromosome has
  //been randomised then reset?
  void testSyncToModelAfterRandomiseReset();
  //14) Does operator== return true for a mutated, reset clone and an original?
  void testChromOperatorEqualsAfterMutateReset();
  //15) Does operator== return true for a randomised, reset clone and an original?
  void testChromOperatorEqualsAfterRandomiseReset();
  //16) Check that chromosomes are still equal after crossover of unmutated clones
  void testNullCrossover();
  //17) Check that conversion from quaternion to Euler angles and back
  //is a null operation
  void testEuler();
  //18) Check that repeated sync to/from model leaves the genotype and phenotype
  //unchanged
  void testRepeatedSync();
  //19) Check that an invalid Compare() returns -1 (i.e. for unequal chromosome lengths)
  void testCompareWithNullChrom();
  //20) Check that Compare() after a single mutation never returns greater than
  //the relative mutation distance (loops over 10000 repeats)
  void testCompareAfterMutate();
  //21) Test operator== with modified _THRESHOLD
  void testOperatorEqualsWithModifiedThreshold();
  //22) Test population constructor
  void testPopulationConstructor();
  //23) Check that population constructor with zero size throws an error
  void testPopulationConstructorZeroSize();
  //24) Check that population constructor with null chromosome throws an error
  void testPopulationConstructorNullChromosome();
  //25) Check that population constructor with null scoring function throws an error
  void testPopulationConstructorNullSF();
  //26) Check that roulette wheel fitness values are calculated OK
  void testPopulationRWFitness();
  //27) Check that roulette wheel selection will select higher scoring genomes on average
  void testPopulationRWSelect();
  //28) Check that population GAstep does not decrease the score of the best genome
  //and that population size remains constant
  void testPopulationGAstep();
  //29) Checks the behaviour of RbtModel::GetChrom if flex data has not been defined
  //Should return a zero length chromosome for rigid model
  void testRbtModelGetChromUndefinedFlexData();
  //30) Checks the behaviour of RbtModel::SetFlexData(NULL)
  //Should return a zero length chromosome for rigid model
  void testRbtModelGetChromNullFlexData();
  //31) Checks that RbtModel::SetFlexData creates the same chromosome
  //as the manual method
  void testRbtModelGetChrom();
  //32) Checks that the variation in randomised values remains in bounds
  //for tethered orientations
  void testRandomiseTetheredOrientation();
  //33) Checks that the variation in randomised values remains in bounds
  //for tethered COM
  void testRandomiseTetheredCOM();
  //34) Checks that the variation in mutated values remains in bounds
  //for tethered orientations
  void testMutateTetheredOrientation();
  //35) Checks that the variation in mutated values remains in bounds
  //for tethered COM
  void testMutateTetheredCOM();
  //36) Checks that the variation in crossed over values remains in bounds
  //for tethered orientations
  void testCrossoverTetheredOrientation();
  //37) Checks that the variation in crossed over values remains in bounds
  //for tethered COM
  void testCrossoverTetheredCOM();
  //38) Checks that the variation in randomised values remains in bounds
  //for tethered dihedral
  void testRandomiseTetheredDihedral();
  //39) Checks that the variation in mutated values remains in bounds
  //for tethered dihedral
  void testMutateTetheredDihedral();
  //40) Checks that the variation in crossed over values remains in bounds
  //for tethered dihedral
  void testCrossoverTetheredDihedral();
  //41) Checks that randomised occupancy is from a rectangular distribution
  void testRandomiseOccupancy();
  //42) Checks that actual occupancy probability matches the desired probability
  void testOccupancyThreshold();
  
private:
  RbtModelPtr m_recep_1koc;
  RbtModelPtr m_lig_1koc;
  RbtDockingSitePtr m_site_1koc;
  RbtChromElementPtr m_chrom_1koc;
  RbtAtomList m_atomList;//All atoms in receptor and ligand
  RbtBiMolWorkSpacePtr m_workSpace;//simple workspace
  RbtSFAggPtr m_SF;//simple scoring function
};
#endif /*RBTCHROMTEST_H_*/
