#include "SearchTest.h"
#include "RbtMdlFileSource.h"
#include "RbtMdlFileSink.h"
#include "RbtPRMFactory.h"
#include "RbtBiMolWorkSpace.h"
#include "RbtVdwIdxSF.h"
#include "RbtVdwIntraSF.h"
#include "RbtCavityGridSF.h"
#include "RbtTransformAgg.h"
#include "RbtGATransform.h"
#include "RbtSimplexTransform.h"
#include "RbtRandPopTransform.h"
#include "RbtSimAnnTransform.h"

CPPUNIT_TEST_SUITE_REGISTRATION( SearchTest );

void SearchTest::setUp() {
    try {
        //Create the docking site, receptor, ligand and solvent objects
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
        //Combine the atom lists of receptor, ligand and solvent
        RbtInt nModels = m_workSpace->GetNumModels();
        for (RbtInt i = 0; i < nModels; i++) {
            RbtAtomList atomList = m_workSpace->GetModel(i)->GetAtomList();
            std::copy(atomList.begin(), atomList.end(),
                    std::back_inserter(m_atomList));
        } 
        //Set up a minimal workspace and scoring function for docking
        m_SF = new RbtSFAgg("SCORE");
        RbtBaseSF* sfInter = new RbtVdwIdxSF("INTER_VDW");
        sfInter->SetParameter(RbtVdwSF::_ECUT, 1.0);
        m_SF->Add(sfInter);
        RbtBaseSF* sfIntra = new RbtVdwIntraSF("INTRA_VDW");
        sfIntra->SetParameter(RbtVdwSF::_ECUT, 1.0);
        m_SF->Add(sfIntra);
        m_workSpace->SetSF(m_SF);
    }
    catch (RbtError& e) {
        cout << e << endl;
    }
}

void SearchTest::tearDown() {
  m_atomList.clear();
  m_SF.SetNull();
  m_workSpace.SetNull();
}

//RMSD calculation between two coordinate lists
RbtDouble SearchTest::rmsd(const RbtCoordList& rc, const RbtCoordList& c) {
  RbtDouble retVal(0.0);
  RbtInt nCoords = rc.size();
  if (c.size() != nCoords) {
    retVal = 999.9;
  }
  else {
    for (RbtInt i = 0; i < nCoords; i++) {
      retVal += Rbt::Length2(rc[i],c[i]);
    }
    retVal = sqrt(retVal/float(nCoords));
  }
  return retVal;
}

void SearchTest::testPRMFactory() {
    CPPUNIT_ASSERT( m_workSpace->GetNumModels() == 6 ); 
}

void SearchTest::testHeavyAtomFactory() {
    RbtAtomRList heavyAtomList;
    //find all the movable heavy atoms in the receptor, ligand and solvent
    if (m_workSpace) {
        RbtCavityGridSF::HeavyAtomFactory atomFactory(m_workSpace->GetModels());
        heavyAtomList = atomFactory.GetAtomList();
    }
    CPPUNIT_ASSERT( heavyAtomList.size() == 42 ); 
}

void SearchTest::testGA() {
    RbtTransformAggPtr spTransformAgg(new RbtTransformAgg());
    RbtBaseTransform* pRandPop = new RbtRandPopTransform();
    RbtBaseTransform* pGA = new RbtGATransform();
    pRandPop->SetParameter(RbtBaseObject::_TRACE, 4);
    pGA->SetParameter(RbtBaseObject::_TRACE, 4);
    spTransformAgg->Add(pRandPop);
    spTransformAgg->Add(pGA);
    m_workSpace->SetTransform(spTransformAgg);
    RbtBool isOK(true);
    try {
        m_workSpace->Run();
    }
    catch (RbtError& e) {
        cout << e.Message() << endl;
        isOK = false;
    }
    CPPUNIT_ASSERT( isOK );    
}

void SearchTest::testSimplex() {
    RbtTransformAggPtr spTransformAgg(new RbtTransformAgg());
    //RbtBaseTransform* pRandPop = new RbtRandPopTransform();
    RbtBaseTransform* pSimplex = new RbtSimplexTransform();
    //pRandPop->SetParameter(RbtBaseObject::_TRACE, 4);
    pSimplex->SetParameter(RbtBaseObject::_TRACE, 1);
    pSimplex->SetParameter(RbtSimplexTransform::_MAX_CALLS, 500);
    pSimplex->SetParameter(RbtSimplexTransform::_NCYCLES, 100);
    pSimplex->SetParameter(RbtSimplexTransform::_STEP_SIZE, 1.0);
    //spTransformAgg->Add(pRandPop);
    spTransformAgg->Add(pSimplex);
    m_workSpace->SetTransform(spTransformAgg);
    RbtBool isOK(true);
    try {
        m_workSpace->Run();
    }
    catch (RbtError& e) {
        cout << e.Message() << endl;
        isOK = false;
    }
    CPPUNIT_ASSERT( isOK );    
}

void SearchTest::testSimAnn() {
    RbtBaseTransform* pSimAnn = new RbtSimAnnTransform();
    pSimAnn->SetParameter(RbtBaseObject::_TRACE, 2);
    pSimAnn->SetParameter(RbtSimAnnTransform::_BLOCK_LENGTH, 100);
    pSimAnn->SetParameter(RbtSimAnnTransform::_NUM_BLOCKS, 50);
    pSimAnn->SetParameter(RbtSimAnnTransform::_START_T, 300.0);
    pSimAnn->SetParameter(RbtSimAnnTransform::_FINAL_T, 50.0);
    pSimAnn->SetParameter(RbtSimAnnTransform::_STEP_SIZE, 0.5);
    pSimAnn->SetParameter(RbtSimAnnTransform::_PARTITION_DIST, 0.0);
    pSimAnn->SetParameter(RbtSimAnnTransform::_PARTITION_FREQ, 0);
    m_workSpace->SetTransform(pSimAnn);
    RbtBool isOK(true);
    try {
        m_workSpace->Run();
    }
    catch (RbtError& e) {
        cout << e.Message() << endl;
        isOK = false;
    }
    delete pSimAnn;
    CPPUNIT_ASSERT( isOK );    
}

void SearchTest::testRestart() {
    RbtTransformAggPtr spTransformAgg(new RbtTransformAgg());
    RbtBaseTransform* pSimplex = new RbtSimplexTransform();
    pSimplex->SetParameter(RbtBaseObject::_TRACE, 1);
    pSimplex->SetParameter(RbtSimplexTransform::_MAX_CALLS, 500);
    pSimplex->SetParameter(RbtSimplexTransform::_NCYCLES, 100);
    pSimplex->SetParameter(RbtSimplexTransform::_STEP_SIZE, 1.0);
    spTransformAgg->Add(pSimplex);
    m_workSpace->SetTransform(spTransformAgg);
    RbtBool isOK(true);
    RbtDouble finalScore(0.0);
    RbtDouble restartScore(0.0);
    try {
        //create an output sink for the minimised ligand/solvent
        RbtMolecularFileSinkPtr spMdlFileSink(new RbtMdlFileSink("restart.sd",RbtModelPtr()));
        m_workSpace->SetSink(spMdlFileSink);
        m_workSpace->Run();
        m_workSpace->Save();
        finalScore = m_workSpace->GetSF()->Score();
        //Reload the receptor, minimised ligand and solvent
        RbtString prmFileName = Rbt::GetRbtFileName("","1YET.prm");
        RbtParameterFileSourcePtr spPrmSource(new RbtParameterFileSource(prmFileName));
        RbtMolecularFileSourcePtr spMdlFileSource(new RbtMdlFileSource("restart.sd",
                    true,true,true));
        //Ligand segment is always called H, solvent will be H2, H3 etc.
        spMdlFileSource->SetSegmentFilterMap
                                  (Rbt::ConvertStringToSegmentMap("H"));
        RbtPRMFactory prmFactory(spPrmSource, m_workSpace->GetDockingSite());
        prmFactory.SetTrace(1);
        m_workSpace->SetReceptor(prmFactory.CreateReceptor());
        m_workSpace->SetLigand(prmFactory.CreateLigand(spMdlFileSource));
        m_workSpace->SetSolvent(prmFactory.CreateSolvent());
        m_workSpace->UpdateModelCoordsFromChromRecords(spMdlFileSource,1);
        restartScore = m_workSpace->GetSF()->Score();
        //RbtStringVariantMap scoreMap;
        //m_workSpace->GetSF()->ScoreMap(scoreMap);
        //for (RbtStringVariantMapConstIter iter = scoreMap.begin(); iter != scoreMap.end(); ++iter) {
        //  cout << iter->first << " = " << iter->second << endl;
        //}
        cout << "Final score = " << finalScore << endl;
        cout << "Restart score = " << restartScore << endl;
    }
    catch (RbtError& e) {
        cout << e.Message() << endl;
        isOK = false;
    }
    CPPUNIT_ASSERT( isOK && (fabs(restartScore - finalScore) < 0.01));    
}


