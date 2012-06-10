//Unit tests for solvent occupancy
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
#ifndef OCCUPANCY_H_
#define OCCUPANCY_H_

#include <cppunit/extensions/HelperMacros.h>

#include "RbtBiMolWorkSpace.h"
#include "RbtModel.h"


class OccupancyTest : public CppUnit::TestFixture  {
    
CPPUNIT_TEST_SUITE( OccupancyTest );
CPPUNIT_TEST( testSolvationSF );
CPPUNIT_TEST( testPolarSF );
CPPUNIT_TEST( testVdwSF );
CPPUNIT_TEST( testVdwSFSolventModes );
CPPUNIT_TEST( testFlexAtomFactoryReceptor );
CPPUNIT_TEST( testFlexAtomFactorySolvent );
CPPUNIT_TEST( testFlexAtomFactoryLigand );
CPPUNIT_TEST_SUITE_END();

public:
  static const RbtDouble TINY = 1E-4;
  //TextFixture methods
  void setUp();
  void tearDown();

  //1) Check that the total desolvation score for receptor / ligand / (disabled) solvent
  //matches the total score with no solvent present
  void testSolvationSF();
  //2) Check that the total polar score for receptor / ligand / (disabled) solvent
  //matches the total score with no solvent present
  void testPolarSF();
  //3) Check that the total vdW score for receptor / ligand / (disabled) solvent
  //matches the total score with no solvent present
  void testVdwSF();
  void testVdwSFSolventModes();
  //5) Checks the results of RbtFlexAtomFactory for receptor flexibility modes
  void testFlexAtomFactoryReceptor();
  //6) Checks the results of RbtFlexAtomFactory for solvent flexibility modes
  void testFlexAtomFactorySolvent();
  //7) Checks the results of RbtFlexAtomFactory for ligand flexibility modes
  void testFlexAtomFactoryLigand();

   
private:
  //Helper functions
  RbtDouble CompareScoresForDisabledAndNoSolvent();
  RbtDouble CompareScoresForSolventModes();
  RbtBool testFlexAtomFactory(RbtModel* pModel, RbtInt expected[9][3]);
  RbtBiMolWorkSpacePtr m_workSpace;
  RbtModelList m_solventList;
};
#endif /*OCCUPANCY_TEST_H_*/
