#include "OccupancyTest.h"
#include "RbtMdlFileSource.h"
#include "RbtMdlFileSink.h"
#include "RbtPRMFactory.h"
#include "RbtSAIdxSF.h"
#include "RbtSetupPolarSF.h"
#include "RbtPolarIdxSF.h"
#include "RbtVdwIdxSF.h"
#include "RbtSFAgg.h"
#include "RbtVariant.h"
#include "RbtFlexAtomFactory.h"
#include "RbtLigandFlexData.h"

CPPUNIT_TEST_SUITE_REGISTRATION( OccupancyTest );

void OccupancyTest::setUp() {
    try {
        //Create the docking site, receptor, and ligand objects
        const RbtString& wsName = "1YET";
        RbtString prmFileName = Rbt::GetRbtFileName("",wsName + ".prm");
        RbtString ligFileName = Rbt::GetRbtFileName("",wsName + "_c.sd");
        RbtString asFileName = Rbt::GetRbtFileName("",wsName + ".as");
        RbtParameterFileSourcePtr spPrmSource(new RbtParameterFileSource(prmFileName));
        RbtMolecularFileSourcePtr spMdlFileSource(new RbtMdlFileSource(ligFileName,
                    true,true,true));
        m_workSpace = new RbtBiMolWorkSpace();
        ifstream istr(asFileName.c_str(),ios_base::in|ios_base::binary);
        m_workSpace->SetDockingSite(new RbtDockingSite(istr));
        istr.close();
        RbtPRMFactory prmFactory(spPrmSource, m_workSpace->GetDockingSite());
        m_workSpace->SetReceptor(prmFactory.CreateReceptor());
        m_workSpace->SetLigand(prmFactory.CreateLigand(spMdlFileSource));
        m_workSpace->SetSolvent(prmFactory.CreateSolvent());
	m_solventList = m_workSpace->GetSolvent();
    }
    catch (RbtError& e) {
        cout << e << endl;
    }
}

void OccupancyTest::tearDown() {
  m_workSpace.SetNull();
  m_solventList.clear();
}


void OccupancyTest::testSolvationSF() {
  RbtBaseSF* sfAgg = new RbtSFAgg("SCORE");
  RbtBaseSF* sfInter = new RbtSFAgg("INTER");
  RbtBaseSF* sfSolv = new RbtSAIdxSF("SOLV");
  sfAgg->Add(sfInter);
  sfInter->Add(sfSolv);
  m_workSpace->SetSF(sfAgg);

  CPPUNIT_ASSERT( CompareScoresForDisabledAndNoSolvent() < TINY);
}

void OccupancyTest::testPolarSF() {
  RbtBaseSF* sfAgg = new RbtSFAgg("SCORE");
  RbtBaseSF* sfInter = new RbtSFAgg("INTER");
  RbtBaseSF* sfSetupPolar = new RbtSetupPolarSF("SETUP");
  RbtBaseSF* sfPolar = new RbtPolarIdxSF("POLAR");
  sfAgg->Add(sfInter);
  //For this term we have to set the interaction range and increment manually
  //These values match those in the standard scoring function
  sfPolar->SetRange(5.31);
  sfPolar->SetParameter(RbtPolarIdxSF::_INCR,3.36);
  sfInter->Add(sfSetupPolar);
  sfInter->Add(sfPolar);
  m_workSpace->SetSF(sfAgg);

  CPPUNIT_ASSERT( CompareScoresForDisabledAndNoSolvent() < TINY);
}

void OccupancyTest::testVdwSF() {
  RbtBaseSF* sfAgg = new RbtSFAgg("SCORE");
  RbtBaseSF* sfInter = new RbtSFAgg("INTER");
  RbtBaseSF* sfVdw = new RbtVdwIdxSF("VDW");
  sfAgg->Add(sfInter);
  sfInter->Add(sfVdw);
  m_workSpace->SetSF(sfAgg);

  CPPUNIT_ASSERT( CompareScoresForDisabledAndNoSolvent() < TINY);
}

void OccupancyTest::testVdwSFSolventModes() {
  RbtBaseSF* sfAgg = new RbtSFAgg("SCORE");
  RbtBaseSF* sfInter = new RbtSFAgg("INTER");
  RbtBaseSF* sfSystem = new RbtSFAgg("SYSTEM");
  RbtBaseSF* sfVdw = new RbtVdwIdxSF("VDW");
  sfVdw->SetTrace(0);
  sfAgg->Add(sfInter);
  sfAgg->Add(sfSystem);
  sfInter->Add(sfVdw);
  m_workSpace->SetSF(sfAgg);

  CPPUNIT_ASSERT( CompareScoresForSolventModes() < TINY);
}

void OccupancyTest::testFlexAtomFactoryReceptor() {
  RbtInt expected[2][3] = { {2034,10,0}, {2044,0,0} };
  RbtModelPtr spReceptor = m_workSpace->GetReceptor();
  RbtBool isTestOK = true;
  //Test two modes:
  //mode 0 is as read from PRM file (flex OH/NH3+)
  //mode 1 is fixed receptor
  for (RbtInt mode = 0; mode < 2; ++mode) {
    if (mode == 1) {
      spReceptor->SetFlexData(NULL);
    }
    RbtFlexAtomFactory flexAtomFactory(spReceptor);
    RbtAtomRList fixedAtomList = flexAtomFactory.GetFixedAtomList();
    RbtAtomRList tetheredAtomList = flexAtomFactory.GetTetheredAtomList();
    RbtAtomRList freeAtomList = flexAtomFactory.GetFreeAtomList();
    RbtInt nFixed = fixedAtomList.size();
    RbtInt nTethered = tetheredAtomList.size();
    RbtInt nFree = freeAtomList.size();
    //cout << "Receptor: #fixed = " << nFixed << "; #tethered = " << nTethered << "; #free = " << nFree << endl;
    //for (RbtAtomRListConstIter iter = tetheredAtomList.begin(); iter != tetheredAtomList.end(); ++iter) {
    //  cout << "Tethered: " << (*iter)->GetAtomName() << ": " << (*iter)->GetUser1Value() << endl;
    //}
    if ( (nFixed != expected[mode][0])
	 || (nTethered != expected[mode][1])
	 || (nFree != expected[mode][2]) ) {
      isTestOK = false;
    }
  }
  CPPUNIT_ASSERT( isTestOK );  
}

void OccupancyTest::testFlexAtomFactorySolvent() {
  //Expected sizes for fixed, tethered and free atom lists
  //for each value of model flexibility mode
  RbtInt expected[9][3] = { {3,0,0}, {0,3,0}, {0,3,0},
			    {0,3,0}, {0,3,0}, {0,3,0},
			    {0,0,3}, {0,0,3}, {0,0,3} };
  RbtBool isTestOK = false;
  RbtModelList solventList = m_workSpace->GetSolvent();
  if (!solventList.empty()) {
    RbtModelPtr solvent = solventList.front();
    isTestOK = testFlexAtomFactory(solventList.front(), expected);
  }
  CPPUNIT_ASSERT( isTestOK );
}

void OccupancyTest::testFlexAtomFactoryLigand() {
  //Expected sizes for fixed, tethered and free atom lists
  //for each value of model flexibility mode
  RbtInt expected[9][3] = { {0,0,44}, {0,0,44}, {0,0,44},
			    {0,0,44}, {0,0,44}, {0,0,44},
			    {0,0,44}, {0,0,44}, {0,0,44} };
  RbtBool isTestOK = false;
  RbtModelPtr spLigand = m_workSpace->GetLigand();
  if (!spLigand.Null()) {
    isTestOK = testFlexAtomFactory(spLigand, expected);
  }
  CPPUNIT_ASSERT( isTestOK );
}

RbtDouble OccupancyTest::CompareScoresForDisabledAndNoSolvent() {
  RbtDouble retVal = 0.0;
  RbtBaseSF* pSF = m_workSpace->GetSF();
  if (pSF) {
    //A) No solvent present
    m_workSpace->RemoveSolvent();
    RbtDouble scoreNoSolvent = pSF->Score();
    //B) With disabled solvent present
    for (RbtModelListIter iter = m_solventList.begin(); iter != m_solventList.end(); ++iter) {
      (*iter)->SetOccupancy(0.0,0.5);
    }
    m_workSpace->SetSolvent(m_solventList);
    RbtDouble scoreDisabledSolvent = pSF->Score();
    //C) With enabled solvent present (not used at present)
    for (RbtModelListIter iter = m_solventList.begin(); iter != m_solventList.end(); ++iter) {
      (*iter)->SetOccupancy(1.0,0.5);
    }
    m_workSpace->SetSolvent(m_solventList);
    RbtDouble scoreEnabledSolvent = pSF->Score();
    retVal = fabs(scoreNoSolvent - scoreDisabledSolvent);
  }
  return retVal;
}

RbtBool OccupancyTest::testFlexAtomFactory(RbtModel* pModel, RbtInt expected[9][3]) {
  RbtBool retVal = true;
  if (pModel) {
    RbtFlexData* pFlexData = pModel->GetFlexData();
    RbtFlexAtomFactory flexAtomFactory;
    for (RbtInt mode = 0; mode < 9; ++mode) {
      RbtChromElement::eMode eTransMode = (RbtChromElement::eMode)(mode / 3);
      RbtChromElement::eMode eRotMode = (RbtChromElement::eMode)(mode % 3);
      RbtString strTransMode = RbtChromElement::ModeToStr(eTransMode);
      RbtString strRotMode = RbtChromElement::ModeToStr(eRotMode);
      pFlexData->SetParameter(RbtLigandFlexData::_TRANS_MODE,strTransMode);
      pFlexData->SetParameter(RbtLigandFlexData::_ROT_MODE,strRotMode);
      pModel->SetFlexData(pFlexData);
      RbtFlexAtomFactory flexAtomFactory(pModel);
      RbtAtomRList fixedAtomList = flexAtomFactory.GetFixedAtomList();
      RbtAtomRList tetheredAtomList = flexAtomFactory.GetTetheredAtomList();
      RbtAtomRList freeAtomList = flexAtomFactory.GetFreeAtomList();
      RbtInt nFixed = fixedAtomList.size();
      RbtInt nTethered = tetheredAtomList.size();
      RbtInt nFree = freeAtomList.size();
      //cout << endl << "Trans = " << strTransMode << ", Rot = " << strRotMode << endl;
      //for (RbtAtomRListConstIter iter = fixedAtomList.begin(); iter != fixedAtomList.end(); ++iter) {
      //  cout << "Fixed: " << (*iter)->GetAtomName() << ": " << (*iter)->GetUser1Value() << endl;
      //}
      //for (RbtAtomRListConstIter iter = tetheredAtomList.begin(); iter != tetheredAtomList.end(); ++iter) {
      //  cout << "Tethered: " << (*iter)->GetAtomName() << ": " << (*iter)->GetUser1Value() << endl;
      //}
      if ( (nFixed != expected[mode][0])
	|| (nTethered != expected[mode][1])
	   || (nFree != expected[mode][2]) ) {
	retVal = false;
      }
    }
  }
  return retVal;
}

RbtDouble OccupancyTest::CompareScoresForSolventModes() {
  RbtDouble retVal = 0.0;
  RbtDouble lastScore = 0.0;
  RbtBaseSF* pSF = m_workSpace->GetSF();
  if (pSF) {
    RbtModelList solventList = m_workSpace->GetSolvent();
    for (RbtInt mode = 0; mode < 9; ++mode) {
      RbtChromElement::eMode eTransMode = (RbtChromElement::eMode)(mode / 3);
      RbtChromElement::eMode eRotMode = (RbtChromElement::eMode)(mode % 3);
      RbtString strTransMode = RbtChromElement::ModeToStr(eTransMode);
      RbtString strRotMode = RbtChromElement::ModeToStr(eRotMode);
      //Change the modes of all except the first solvent model, to ensure a mix of
      //tethered / free solvent
      for (RbtModelListIter iter = solventList.begin()+1; iter != solventList.end(); ++iter) {
	RbtFlexData* pFlexData = (*iter)->GetFlexData();
	pFlexData->SetParameter(RbtLigandFlexData::_TRANS_MODE,strTransMode);
	pFlexData->SetParameter(RbtLigandFlexData::_ROT_MODE,strRotMode);
	(*iter)->SetFlexData(pFlexData);
      }
      m_workSpace->RemoveSolvent();
      m_workSpace->SetSolvent(solventList);
      RbtDouble score = pSF->Score();
      //cout << "mode = " << mode << "; score = " << score << endl;
      if (mode > 0) {
	RbtDouble diff = fabs(score-lastScore);
	retVal = std::max(diff, retVal);
      }
      lastScore = score;
      //RbtStringVariantMap scoreMap;
      //pSF->ScoreMap(scoreMap);
      //for (RbtStringVariantMapConstIter vIter = scoreMap.begin(); vIter != scoreMap.end(); vIter++) {
      //cout << (*vIter).first << " = " << (*vIter).second << endl;
      //}
    }
  }
  return retVal;
}
