//Unit tests for rDock search algorithms
//
//Required input files:
//1YET.prm  rDock receptor file
//1YET.psf  Receptor topology file
//1YET.crd  Receptor coordinate file
//1YET_c.sd Ligand coordinate file
//1YET.as   Docking site
//
//Required environment:
//Make sure the above files are colocated in a single directory
//and define RBT_HOME env. variable to point at this directory
#ifndef SEARCHTEST_H_
#define SEARCHTEST_H_

#include <cppunit/extensions/HelperMacros.h>

#include "RbtChrom.h"
#include "RbtModel.h"
#include "RbtDockingSite.h"
#include "RbtBiMolWorkSpace.h"
#include "RbtSFAgg.h"

class SearchTest : public CppUnit::TestFixture  {
    
CPPUNIT_TEST_SUITE( SearchTest );
CPPUNIT_TEST( testPRMFactory );
CPPUNIT_TEST( testHeavyAtomFactory );
CPPUNIT_TEST( testGA );
CPPUNIT_TEST( testSimplex );
CPPUNIT_TEST( testSimAnn );
CPPUNIT_TEST( testRestart );
CPPUNIT_TEST_SUITE_END();

public:
  static const RbtDouble TINY = 1E-4;
  //TextFixture methods
  void setUp();
  void tearDown();

  //rdock helper methods
  //RMSD calculation between two coordinate lists
  RbtDouble rmsd(const RbtCoordList& rc, const RbtCoordList& c);

  //1 Check that receptor, ligand and solvent models are loaded into workspace
  //Should be 6 models in total (4 solvent)
  void testPRMFactory();
  //2 Check RbtFlexDataVisitor subclass correctly identifies movable heavy atoms in cavity
  void testHeavyAtomFactory();
  //3 Run a sample GA
  void testGA();
  //4 Run a sample Simplex
  void testSimplex();
  //5 Run a sample simulated annealing
  void testSimAnn();
  //6 Check we can reload solvent coords from ligand SD file
  void testRestart();
   
private:
  RbtAtomList m_atomList;//All atoms in receptor, ligand and solvent
  RbtBiMolWorkSpacePtr m_workSpace;//simple workspace
  RbtSFAggPtr m_SF;//simple scoring function
};
#endif /*SEARCHTEST_H_*/
