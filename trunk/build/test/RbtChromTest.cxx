    #include "RbtChromTest.h"
    #include "RbtPRMFactory.h"
    #include "RbtLigandFlexData.h"
    #include "RbtReceptorFlexData.h"
    #include "RbtChromFactory.h"
    #include "RbtChromPositionElement.h"
    #include "RbtChromDihedralElement.h"
    #include "RbtChromOccupancyElement.h"
    #include "RbtMdlFileSource.h"
    #include "RbtEuler.h"
    #include "RbtRand.h"
    #include "RbtSFAgg.h"
    #include "RbtVdwIdxSF.h"
    #include "RbtVdwIntraSF.h"
    #include "RbtBiMolWorkSpace.h"
    
    CPPUNIT_TEST_SUITE_REGISTRATION( RbtChromTest );
    
    void RbtChromTest::setUp() {
        try {
            //Create a receptor model, ligand model, and docking site
            const RbtString& wsName = "1koc";
            RbtString prmFileName = Rbt::GetRbtFileName("data/receptors",wsName + ".prm");
            RbtString ligFileName = Rbt::GetRbtFileName("",wsName + "_c.sd");
            RbtString asFileName = Rbt::GetRbtFileName("data/grids",wsName + ".as");
            RbtParameterFileSourcePtr spPrmSource(new RbtParameterFileSource(prmFileName));
            RbtMolecularFileSourcePtr spMdlFileSource(new RbtMdlFileSource(ligFileName,
                        true,true,true));
            ifstream istr(asFileName.c_str(),ios_base::in|ios_base::binary);
            m_site_1koc = new RbtDockingSite(istr);
            istr.close();
            RbtPRMFactory prmFactory(spPrmSource);
            m_recep_1koc = prmFactory.CreateReceptor();
            m_lig_1koc = prmFactory.CreateLigand(spMdlFileSource);
            //Create a ligand chromosome for free docking
            m_chrom_1koc = new RbtChrom();
            RbtFlexDataPtr ligFlexData
                    = new RbtLigandFlexData(m_site_1koc);
            ligFlexData->SetModel(m_lig_1koc);
            RbtChromFactory chromFactory1;
            ligFlexData->Accept(chromFactory1);
            m_chrom_1koc->Add(chromFactory1.GetChrom());
            //Create a receptor chromosome for flexible terminal bonds
            RbtFlexDataPtr recFlexData
                    = new RbtReceptorFlexData(m_site_1koc);
            recFlexData->SetModel(m_recep_1koc);
            RbtChromFactory chromFactory2;
            recFlexData->Accept(chromFactory2);
            m_chrom_1koc->Add(chromFactory2.GetChrom());
            //Combine the atom lists of ligand and receptor
            m_atomList = m_lig_1koc->GetAtomList();
            RbtAtomList recepAtomList = m_recep_1koc->GetAtomList();
            std::copy(recepAtomList.begin(), recepAtomList.end(),
                        std::back_inserter(m_atomList));
        }
        catch (RbtError& e) {
            cout << e << endl;
        }
    }
    
    void RbtChromTest::tearDown() {
      m_recep_1koc.SetNull();
      m_lig_1koc.SetNull();
      m_site_1koc.SetNull();
      m_chrom_1koc.SetNull();
      m_atomList.clear();
      m_SF.SetNull();
      m_workSpace.SetNull();
    }
    
    //extra setup/teardown for building a simple scoring function and workspace
    void RbtChromTest::setupWorkSpace() {
        //Set up a minimal workspace and scoring function for docking
        m_SF = new RbtSFAgg("SCORE");
        RbtBaseSF* sfInter = new RbtVdwIdxSF("INTER_VDW");
        sfInter->SetParameter(RbtVdwSF::_ECUT, 1.0);
        m_SF->Add(sfInter);
        RbtBaseSF* sfIntra = new RbtVdwIntraSF("INTRA_VDW");
        sfIntra->SetParameter(RbtVdwSF::_ECUT, 1.0);
        m_SF->Add(sfIntra);
        m_workSpace = new RbtBiMolWorkSpace();
        m_workSpace->SetSF(m_SF);
        m_workSpace->SetDockingSite(m_site_1koc);
        m_workSpace->SetReceptor(m_recep_1koc);
        m_workSpace->SetLigand(m_lig_1koc);    
    }
    
    //RMSD calculation between two coordinate lists
    RbtDouble RbtChromTest::rmsd(const RbtCoordList& rc, const RbtCoordList& c) {
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
    
    
    void RbtChromTest::testChromFactory() {
        CPPUNIT_ASSERT(m_chrom_1koc.Ptr() != NULL);    
    }
    
    void RbtChromTest::testCloneDestructor() {
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        CPPUNIT_ASSERT(m_chrom_1koc.Ptr() != NULL);    
    }
    
    void RbtChromTest::testCloneLengthEquals() {
       RbtChromElementPtr clone = m_chrom_1koc->clone();
       CPPUNIT_ASSERT(m_chrom_1koc->GetLength() == clone->GetLength());
    }
    
    void RbtChromTest::testChromOperatorEquals() {
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        CPPUNIT_ASSERT(*m_chrom_1koc == *clone);
    }
    
    void RbtChromTest::testChromOperatorEqualsAfterMutate() {
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        clone->Mutate(1.0);
        CPPUNIT_ASSERT(*m_chrom_1koc != *clone);
    }
    
    void RbtChromTest::testSyncToModel() {
        RbtCoordList coordsBefore, coordsAfter;
        Rbt::GetCoordList(m_atomList, coordsBefore);
        m_chrom_1koc->SyncToModel();
        Rbt::GetCoordList(m_atomList, coordsAfter);
        RbtDouble rms = rmsd(coordsBefore, coordsAfter);
        CPPUNIT_ASSERT(rms < TINY);     
    }
    
    void RbtChromTest::testSyncToModelAfterMutate() {
        RbtCoordList coordsBefore, coordsAfter;
        Rbt::GetCoordList(m_atomList, coordsBefore);
        m_chrom_1koc->Mutate(1.0);
        m_chrom_1koc->SyncToModel();
        Rbt::GetCoordList(m_atomList, coordsAfter);
        RbtDouble rms = rmsd(coordsBefore, coordsAfter);
        CPPUNIT_ASSERT(rms > TINY);     
    }
    
    void RbtChromTest::testSyncToModelAfterReset() {
        RbtCoordList coordsBefore, coordsAfter;
        Rbt::GetCoordList(m_atomList, coordsBefore);
        m_chrom_1koc->Mutate(1.0);
        m_chrom_1koc->SyncToModel();
        m_chrom_1koc->Reset();
        m_chrom_1koc->SyncToModel();
        Rbt::GetCoordList(m_atomList, coordsAfter);
        RbtDouble rms = rmsd(coordsBefore, coordsAfter);
        CPPUNIT_ASSERT(rms < TINY);     
    }
    
    void RbtChromTest::testCloneReset() {
        RbtCoordList coordsBefore, coordsAfter;
        Rbt::GetCoordList(m_atomList, coordsBefore);
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        clone->Mutate(1.0);
        clone->SyncToModel();
        clone->Reset();
        clone->SyncToModel();
        Rbt::GetCoordList(m_atomList, coordsAfter);
        RbtDouble rms = rmsd(coordsBefore, coordsAfter);
        CPPUNIT_ASSERT(rms < TINY);
    }
    
    void RbtChromTest::testMutatedCloneReset() {
        RbtCoordList coordsBefore, coordsAfter;
        Rbt::GetCoordList(m_atomList, coordsBefore);
        //Mutate the original before creating the clone
        m_chrom_1koc->Mutate(1.0);
        m_chrom_1koc->SyncToModel();
        RbtChromElement* clone = m_chrom_1koc->clone();
        clone->Mutate(1.0);
        clone->SyncToModel();
        clone->Reset();
        clone->SyncToModel();
        Rbt::GetCoordList(m_atomList, coordsAfter);
        RbtDouble rms = rmsd(coordsBefore, coordsAfter);
        CPPUNIT_ASSERT(rms < TINY);
        delete clone;
    }
    
    void RbtChromTest::testChromOperatorEqualsAfterRandomise() {
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        clone->Randomise();
        CPPUNIT_ASSERT(*m_chrom_1koc != *clone);
    }
    
    void RbtChromTest::testSyncToModelAfterRandomise() {
        RbtCoordList coordsBefore, coordsAfter;
        Rbt::GetCoordList(m_atomList, coordsBefore);
        m_chrom_1koc->Randomise();
        m_chrom_1koc->SyncToModel();
        Rbt::GetCoordList(m_atomList, coordsAfter);
        RbtDouble rms = rmsd(coordsBefore, coordsAfter);
        CPPUNIT_ASSERT(rms > TINY);
    }
    
    void RbtChromTest::testSyncToModelAfterRandomiseReset() {
        RbtCoordList coordsBefore, coordsAfter;
        Rbt::GetCoordList(m_atomList, coordsBefore);
        m_chrom_1koc->Randomise();
        m_chrom_1koc->SyncToModel();
        m_chrom_1koc->Reset();
        m_chrom_1koc->SyncToModel();
        Rbt::GetCoordList(m_atomList, coordsAfter);
        RbtDouble rms = rmsd(coordsBefore, coordsAfter);
        CPPUNIT_ASSERT(rms < TINY);         
    }
    
    void RbtChromTest::testChromOperatorEqualsAfterMutateReset() {
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        clone->Mutate(1.0);
        clone->Reset();
        CPPUNIT_ASSERT(*m_chrom_1koc == *clone);
    }
    
    void RbtChromTest::testChromOperatorEqualsAfterRandomiseReset() {
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        clone->Randomise();
        clone->Reset();
        CPPUNIT_ASSERT(*m_chrom_1koc == *clone);
    }
    
    void RbtChromTest::testNullCrossover() {
        RbtChromElementPtr clone1 = m_chrom_1koc->clone();
        RbtChromElementPtr clone2 = m_chrom_1koc->clone();
        RbtChromElementPtr clone3 = m_chrom_1koc->clone();
        RbtChromElementPtr clone4 = m_chrom_1koc->clone();
        RbtBool isEqualBeforeXOver = (*clone1 == *clone2)
                                && (*clone1 == *clone3)
                                && (*clone1 == *clone4);
        Rbt::Crossover(clone1, clone2, clone3, clone4);
        RbtBool isEqualAfterXOver = (*clone1 == *clone2)
                                && (*clone1 == *clone3)
                                && (*clone1 == *clone4);
        CPPUNIT_ASSERT(isEqualBeforeXOver && isEqualAfterXOver);
    }
    
    void RbtChromTest::testEuler() {
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        RbtAtomList atomList = m_lig_1koc->GetAtomList();
        RbtPrincipalAxes prAxes = Rbt::GetPrincipalAxes(atomList);
        RbtPrincipalAxes cartesianAxes;
        RbtQuat q = Rbt::GetQuatFromAlignAxes(cartesianAxes, prAxes);
        RbtEuler euler(q);
        RbtQuat q1 = euler.ToQuat();
        CPPUNIT_ASSERT( Rbt::Length(q1-q) < TINY );
    }
    
    void RbtChromTest::testRepeatedSync() {
        //Keep clone unchanged as a reference
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        RbtCoordList coordsBefore, coordsAfter;
        Rbt::GetCoordList(m_atomList, coordsBefore);
        for (RbtInt i = 0; i < 100; i++) {
            m_chrom_1koc->SyncToModel();
            m_chrom_1koc->SyncFromModel();
        }
        Rbt::GetCoordList(m_atomList, coordsAfter);
        RbtDouble rms = rmsd(coordsBefore, coordsAfter);
        RbtDouble cmp = m_chrom_1koc->Compare(*clone);
        CPPUNIT_ASSERT( (rms < TINY) && (cmp < TINY) );
    }
    
    void RbtChromTest::testCompareWithNullChrom() {
        RbtChromElementPtr nullChrom = new RbtChrom();
        RbtDouble cmp = m_chrom_1koc->Compare(*nullChrom);
        CPPUNIT_ASSERT( (cmp < 0.0) && (*m_chrom_1koc != *nullChrom) );
    }
    
    void RbtChromTest::testCompareAfterMutate() {
        RbtBool isOK(true);
        RbtRand& rand = Rbt::GetRbtRand();
        for (RbtInt i = 0; (i < 10000) && isOK; i++) {
            RbtChromElementPtr clone = m_chrom_1koc->clone();
            RbtDouble mutationDistance = rand.GetRandom01();
            clone->Mutate(mutationDistance);
            RbtDouble cmp = m_chrom_1koc->Compare(*clone);
            isOK = (cmp >= 0.0) && (cmp <= mutationDistance);
        }
        CPPUNIT_ASSERT( isOK );
    }
    
    void RbtChromTest::testOperatorEqualsWithModifiedThreshold() {
        RbtChromElementPtr clone = m_chrom_1koc->clone();
        clone->SyncToModel();
        clone->SyncFromModel();
        //Very fine threshold - should fail
        RbtChromElement::_THRESHOLD = 1E-20;
        RbtBool isEqual20 = (*m_chrom_1koc == *clone);
        //Fine threshold, but should pass
        RbtChromElement::_THRESHOLD = 1E-10;
        RbtBool isEqual10 = (*m_chrom_1koc == *clone);
        //Crude threshold - should pass even after a single mutation
        RbtChromElement::_THRESHOLD = 1.0;
        clone->Mutate(0.5);
        RbtBool isEqualAfterMutate1 = (*m_chrom_1koc == *clone);
        CPPUNIT_ASSERT( !isEqual20 && isEqual10 && isEqualAfterMutate1 );
    }
    
    void RbtChromTest::testPopulationConstructor() {
        setupWorkSpace();
        RbtInt popSize = 100;
        RbtPopulationPtr pop = new RbtPopulation(m_chrom_1koc, popSize, m_SF);
        CPPUNIT_ASSERT( (pop->GetMaxSize() == popSize)
                        && (pop->GetActualSize() == popSize) );    
    }
    
    void RbtChromTest::testPopulationConstructorZeroSize() {
        setupWorkSpace();
        RbtInt popSize = 0;
        RbtPopulationPtr pop = new RbtPopulation(m_chrom_1koc, popSize, m_SF);
    }
    
    void RbtChromTest::testPopulationConstructorNullChromosome() {
        setupWorkSpace();
        RbtInt popSize = 100;
        RbtPopulationPtr pop = new RbtPopulation(NULL, popSize, m_SF);
    }
    
    void RbtChromTest::testPopulationConstructorNullSF() {
        setupWorkSpace();
        RbtInt popSize = 100;
        RbtPopulationPtr pop = new RbtPopulation(m_chrom_1koc, popSize, NULL);
    }
    
    void RbtChromTest::testPopulationRWFitness() {
        setupWorkSpace();
        RbtInt popSize = 100;
        RbtPopulationPtr pop = new RbtPopulation(m_chrom_1koc, popSize, m_SF);
        const RbtGenomeList& genomeList = pop->GetGenomeList();
        RbtDouble lastValue = 0.0;
        RbtBool isAscending = true;
        for (RbtGenomeListConstIter iter = genomeList.begin();
                    (iter != genomeList.end()) && isAscending; ++iter) {
            RbtDouble value = (*iter)->GetRWFitness();
            //cout << (*iter)->GetScore() << "\t" << (*iter)->GetRWFitness() << endl;
            isAscending = (value >= lastValue);
            lastValue = value;
        }
        CPPUNIT_ASSERT( isAscending && (fabs(lastValue-1.0) < TINY));
    }
    
    void RbtChromTest::testPopulationRWSelect() {
        setupWorkSpace();
        RbtInt popSize = 100;
        RbtPopulationPtr pop = new RbtPopulation(m_chrom_1koc, popSize, m_SF);
        RbtDouble popMean = pop->GetScoreMean();
        RbtDouble selectionMean(0.0);
        RbtInt nToSelect(100);
        for (RbtInt i = 0; i < nToSelect; i++) {
            RbtGenomePtr genome = pop->RouletteWheelSelect();
            RbtDouble score = genome->GetScore();
            selectionMean += score;
            //cout << "RWSelect score=" << score << endl;
        }
        selectionMean /= nToSelect;
        //cout << "popMean = " << popMean << "; selectionMean = " << selectionMean << endl;
        CPPUNIT_ASSERT( selectionMean > popMean );
    }
    
    void RbtChromTest::testPopulationGAstep() {
        setupWorkSpace();
        RbtInt popSize = 100;
        RbtInt nReplicates = 50;
        RbtInt nIter = 100;
        RbtDouble equalityThreshold = 1.0E-2;
        RbtDouble relStepSize = 1.0;
        RbtDouble pcross = 0.4;
        RbtBool xovermut = true;
        RbtBool cmutate = false;
        RbtPopulationPtr pop = new RbtPopulation(m_chrom_1koc, popSize, m_SF);
        RbtDouble lastScore = pop->Best()->GetScore();
        RbtBool isOK = true;
        try {
            for (RbtInt i = 0; (i < nIter) && isOK; ++i) {
                pop->GAstep(nReplicates, relStepSize, equalityThreshold, pcross,
                            xovermut, cmutate);
                RbtDouble score = pop->Best()->GetScore();
                isOK = (score >= lastScore);
                lastScore = score;
            }
        }
        catch (RbtError& e) {
            cout << e.Message() << endl;
            isOK = false;
        }
        CPPUNIT_ASSERT( isOK );    
    }
    
    void RbtChromTest::testRbtModelGetChromUndefinedFlexData() {
        RbtChromElementPtr chrom = new RbtChrom();
        chrom->Add(m_lig_1koc->GetChrom());
        chrom->Add(m_recep_1koc->GetChrom());
        CPPUNIT_ASSERT( chrom->GetLength() == 0);
    }
    
    void RbtChromTest::testRbtModelGetChromNullFlexData() {
        m_lig_1koc->SetFlexData(NULL);
        m_recep_1koc->SetFlexData(NULL);
        RbtChromElementPtr chrom = new RbtChrom();
        chrom->Add(m_lig_1koc->GetChrom());
        chrom->Add(m_recep_1koc->GetChrom());
        CPPUNIT_ASSERT( chrom->GetLength() == 0);
    }
    
    void RbtChromTest::testRbtModelGetChrom() {
        m_lig_1koc->SetFlexData(
                        new RbtLigandFlexData(m_site_1koc));
        m_recep_1koc->SetFlexData(
                        new RbtReceptorFlexData(m_site_1koc));
        RbtChromElementPtr chrom = new RbtChrom();
        chrom->Add(m_lig_1koc->GetChrom());
        chrom->Add(m_recep_1koc->GetChrom());
        CPPUNIT_ASSERT( *chrom == *m_chrom_1koc);
    }
    
    void RbtChromTest::testRandomiseTetheredOrientation() {
        RbtDouble transStepSize = 0.1;
        RbtDouble rotStepSize = 10.0;
        RbtChromElement::eMode transMode = RbtChromElement::FIXED;
        RbtChromElement::eMode rotMode = RbtChromElement::TETHERED;
        RbtDouble maxTrans = 1.0;
        RbtDouble maxRot = 45.0;
        rotStepSize *= M_PI/180.0;
        maxRot *= M_PI/180.0;
        RbtChromElementPtr chrom = new RbtChromPositionElement(
                                        m_lig_1koc,
                                        m_site_1koc,
                                        transStepSize,
                                        rotStepSize,
                                        transMode,
                                        rotMode,
                                        maxTrans,
                                        maxRot);
        RbtDouble meanDiff, minDiff, maxDiff;
        //Measure 1M randomisations of the tethered orientation
        measureRandOrMutateDiff(chrom, 1000000, false,
                                meanDiff, minDiff, maxDiff);
        //The expected mean rotation should be half the maximum step size
        //as CompareVector always returns a positive difference.
        RbtDouble expectedMean = 0.5 * maxRot / rotStepSize;
        //Check that min, max rotation are in range [0,maxRot]
        //and that mean rotation is within 1% of the expected mean.
        CPPUNIT_ASSERT ( (minDiff >= 0.0)
                      && (maxDiff <= (2.0 * expectedMean))
                      && ((meanDiff / expectedMean) > 0.99)
                      && ((meanDiff / expectedMean) < 1.01));
    }
    
    void RbtChromTest::testRandomiseTetheredCOM() {
        RbtDouble transStepSize = 0.1;
        RbtDouble rotStepSize = 10.0;
        RbtChromElement::eMode transMode = RbtChromElement::TETHERED;
        RbtChromElement::eMode rotMode = RbtChromElement::FIXED;
        RbtDouble maxTrans = 1.0;
        RbtDouble maxRot = 45.0;
        rotStepSize *= M_PI/180.0;
        maxRot *= M_PI/180.0;
        RbtChromElementPtr chrom = new RbtChromPositionElement(
                                        m_lig_1koc,
                                        m_site_1koc,
                                        transStepSize,
                                        rotStepSize,
                                        transMode,
                                        rotMode,
                                        maxTrans,
                                        maxRot);
        RbtDouble meanDiff, minDiff, maxDiff;
        //Measure 1M randomisations of the tethered COM
        measureRandOrMutateDiff(chrom, 1000000, false,
                                meanDiff, minDiff, maxDiff);
        //The expected mean displacement should be half the maximum step size
        //as CompareVector always returns a positive difference.
        RbtDouble expectedMean = 0.5 * maxTrans / transStepSize;
        //Check that min, max displacement are in range [0,maxTrans]
        //and that mean rotation is within 1% of the expected mean.
        CPPUNIT_ASSERT ( (minDiff >= 0.0)
                      && (maxDiff <= (2.0 * expectedMean))
                      && ((meanDiff / expectedMean) > 0.99)
                      && ((meanDiff / expectedMean) < 1.01));
    }
    
    void RbtChromTest::testRandomiseTetheredDihedral() {
        RbtDouble stepSize = 10.0;
        RbtChromElement::eMode mode = RbtChromElement::TETHERED;
        RbtDouble maxDelta = 45.0;
        RbtBondList rotBondList = Rbt::GetBondList(m_lig_1koc->GetBondList(), Rbt::isBondRotatable());
        RbtAtomList noTetheredAtoms;
        RbtBondListConstIter iter = rotBondList.begin();
        RbtChromElementPtr chrom = new RbtChromDihedralElement(*iter,
                                            noTetheredAtoms,
                                            stepSize,
                                            mode,
                                            maxDelta);
        RbtDouble meanDiff, minDiff, maxDiff;
        //Measure 1M randomisations of the tethered dihedrals
        measureRandOrMutateDiff(chrom, 1000000, false,
                                meanDiff, minDiff, maxDiff);
        //The expected mean displacement should be half the maximum step size
        //as CompareVector always returns a positive difference.
        RbtDouble expectedMean = 0.5 * maxDelta / stepSize;
        //Check that min, max displacement are in range [0,maxTrans]
        //and that mean rotation is within 1% of the expected mean.
        CPPUNIT_ASSERT ( (minDiff >= 0.0)
                      && (maxDiff <= (2.0 * expectedMean))
                      && ((meanDiff / expectedMean) > 0.99)
                      && ((meanDiff / expectedMean) < 1.01));
    }
    
    void RbtChromTest::testMutateTetheredOrientation() {
        RbtDouble transStepSize = 0.1;
        RbtDouble rotStepSize = 10.0;
        RbtChromElement::eMode transMode = RbtChromElement::FIXED;
        RbtChromElement::eMode rotMode = RbtChromElement::TETHERED;
        RbtDouble maxTrans = 1.0;
        RbtDouble maxRot = 45.0;
        rotStepSize *= M_PI/180.0;
        maxRot *= M_PI/180.0;
        RbtChromElementPtr chrom = new RbtChromPositionElement(
                                        m_lig_1koc,
                                        m_site_1koc,
                                        transStepSize,
                                        rotStepSize,
                                        transMode,
                                        rotMode,
                                        maxTrans,
                                        maxRot);
        RbtDouble meanDiff, minDiff, maxDiff;
        measureRandOrMutateDiff(chrom, 1000000, true,
                                meanDiff, minDiff, maxDiff);
        CPPUNIT_ASSERT ( (minDiff >= 0.0)
                      && (maxDiff <= 1.01 * (maxRot / rotStepSize)));
    }
    
    void RbtChromTest::testMutateTetheredCOM() {
        RbtDouble transStepSize = 0.1;
        RbtDouble rotStepSize = 10.0;
        RbtChromElement::eMode transMode = RbtChromElement::TETHERED;
        RbtChromElement::eMode rotMode = RbtChromElement::FIXED;
        RbtDouble maxTrans = 1.0;
        RbtDouble maxRot = 45.0;
        rotStepSize *= M_PI/180.0;
        maxRot *= M_PI/180.0;
        RbtChromElementPtr chrom = new RbtChromPositionElement(
                                        m_lig_1koc,
                                        m_site_1koc,
                                        transStepSize,
                                        rotStepSize,
                                        transMode,
                                        rotMode,
                                        maxTrans,
                                        maxRot);
        RbtDouble meanDiff, minDiff, maxDiff;
        measureRandOrMutateDiff(chrom, 1000000, true,
                                meanDiff, minDiff, maxDiff);
        CPPUNIT_ASSERT ( (minDiff >= 0.0)
                      && (maxDiff <= 1.01 * (maxTrans / transStepSize)));
    }
    
    void RbtChromTest::testMutateTetheredDihedral() {
        RbtDouble stepSize = 10.0;
        RbtChromElement::eMode mode = RbtChromElement::TETHERED;
        RbtDouble maxDelta = 45.0;
        RbtBondList rotBondList = Rbt::GetBondList(m_lig_1koc->GetBondList(), Rbt::isBondRotatable());
        RbtAtomList noTetheredAtoms;
        RbtBondListConstIter iter = rotBondList.begin();
        RbtChromElementPtr chrom = new RbtChromDihedralElement(*iter,
                                            noTetheredAtoms,
                                            stepSize,
                                            mode,
                                            maxDelta);
        RbtDouble meanDiff, minDiff, maxDiff;
        measureRandOrMutateDiff(chrom, 1000000, true,
                                meanDiff, minDiff, maxDiff);
        CPPUNIT_ASSERT ( (minDiff >= 0.0)
                      && (maxDiff <= 1.01 * (maxDelta / stepSize)));
    }
    
    void RbtChromTest::testCrossoverTetheredOrientation() {
        RbtDouble transStepSize = 0.1;
        RbtDouble rotStepSize = 10.0;
        RbtChromElement::eMode transMode = RbtChromElement::FIXED;
        RbtChromElement::eMode rotMode = RbtChromElement::TETHERED;
        RbtDouble maxTrans = 1.0;
        RbtDouble maxRot = 45.0;
        rotStepSize *= M_PI/180.0;
        maxRot *= M_PI/180.0;
        RbtChromElementPtr chrom = new RbtChromPositionElement(
                                        m_lig_1koc,
                                        m_site_1koc,
                                        transStepSize,
                                        rotStepSize,
                                        transMode,
                                        rotMode,
                                        maxTrans,
                                        maxRot);
        RbtDouble meanDiff, minDiff, maxDiff;
        measureCrossoverDiff(chrom, 100, meanDiff, minDiff, maxDiff);
        CPPUNIT_ASSERT ( (minDiff >= 0.0)
                      && (maxDiff <= 1.01 * (maxRot / rotStepSize)));
    }
    
    void RbtChromTest::testCrossoverTetheredCOM() {
        RbtDouble transStepSize = 0.1;
        RbtDouble rotStepSize = 10.0;
        RbtChromElement::eMode transMode = RbtChromElement::TETHERED;
        RbtChromElement::eMode rotMode = RbtChromElement::FIXED;
        RbtDouble maxTrans = 1.0;
        RbtDouble maxRot = 45.0;
        rotStepSize *= M_PI/180.0;
        maxRot *= M_PI/180.0;
        RbtChromElementPtr chrom = new RbtChromPositionElement(
                                        m_lig_1koc,
                                        m_site_1koc,
                                        transStepSize,
                                        rotStepSize,
                                        transMode,
                                        rotMode,
                                        maxTrans,
                                        maxRot);
        RbtDouble meanDiff, minDiff, maxDiff;
        measureCrossoverDiff(chrom, 100, meanDiff, minDiff, maxDiff);
        CPPUNIT_ASSERT ( (minDiff >= 0.0)
                      && (maxDiff <= 1.01 * (maxTrans / transStepSize)));
    }
    
    void RbtChromTest::testCrossoverTetheredDihedral() {
        RbtDouble stepSize = 10.0;
        RbtChromElement::eMode mode = RbtChromElement::TETHERED;
        RbtDouble maxDelta = 45.0;
        RbtBondList rotBondList = Rbt::GetBondList(m_lig_1koc->GetBondList(), Rbt::isBondRotatable());
        RbtAtomList noTetheredAtoms;
        RbtBondListConstIter iter = rotBondList.begin();
        RbtChromElementPtr chrom = new RbtChromDihedralElement(*iter,
                                            noTetheredAtoms,
                                            stepSize,
                                            mode,
                                            maxDelta);
        RbtDouble meanDiff, minDiff, maxDiff;
        measureCrossoverDiff(chrom, 100, meanDiff, minDiff, maxDiff);
        CPPUNIT_ASSERT ( (minDiff >= 0.0)
                      && (maxDiff <= 1.01 * (maxDelta / stepSize)));
    }

    void RbtChromTest::testRandomiseOccupancy() {
	RbtDouble stepSize = 0.1;
	RbtDouble threshold = 0.5;
	RbtChromElementPtr chrom = new RbtChromOccupancyElement(
					m_lig_1koc,
					stepSize,
					threshold);
	RbtDouble meanDiff, minDiff, maxDiff;
        measureRandOrMutateDiff(chrom, 1000000, false,
				meanDiff, minDiff, maxDiff);
	//Randomised occupancy should be a rectangular distribution
	//between 0 and 1, mean = 0.5.
	minDiff *= stepSize;
	meanDiff *= stepSize;
	maxDiff *= stepSize;
        CPPUNIT_ASSERT (       (fabs(minDiff) < 0.01)
			    && (fabs(meanDiff - 0.5) < 0.01)
			    && (fabs(maxDiff  - 1.0) < 0.01));
    }

    void RbtChromTest::testOccupancyThreshold() {
	RbtDouble stepSize = 0.1;
	RbtDouble occupancyProb = 0.7;
	RbtDouble threshold = 1.0 - occupancyProb;
	RbtChromElementPtr chrom = new RbtChromOccupancyElement(
					m_lig_1koc,
					stepSize,
					threshold);
	RbtInt nTrials = 1000000;
	RbtInt nEnabled = 0;
	RbtAtomList atomList = m_lig_1koc->GetAtomList();
	RbtAtom* pAtom0 = atomList.front();
	//Measure whether the actual probability of enabling the model atoms
	//is equal to the desired probability
	for (RbtInt i = 0; i < nTrials; i++) {
	  chrom->Randomise();
	  chrom->SyncToModel();
	  if (pAtom0->GetEnabled()) {
	    nEnabled++;
	  }
	}
	RbtDouble enabledProb = (RbtDouble)nEnabled / (RbtDouble)nTrials;
        CPPUNIT_ASSERT ( (fabs(enabledProb - occupancyProb) < 0.01));
    }



    void RbtChromTest::measureRandOrMutateDiff(RbtChromElement* chrom,
                                            RbtInt nTrials,
                                            RbtBool bMutate,
                                            RbtDouble& meanDiff,
                                            RbtDouble& minDiff,
                                            RbtDouble& maxDiff) {
        RbtDoubleList refVec;
        chrom->GetVector(refVec);
        meanDiff = 0.0;
        minDiff = 0.0;
        maxDiff = 0.0;
        for (RbtInt iTrial = 0; iTrial < nTrials; iTrial++) {
            if (bMutate) {
                chrom->Mutate(1.0);
            }
            else {
                chrom->Randomise();
            }
            RbtInt i(0);
            RbtDouble diff = chrom->CompareVector(refVec,i);
            meanDiff += diff;
            if (iTrial == 0) {
                minDiff = diff;
                maxDiff = diff;
            }
            else {
                minDiff = std::min(minDiff, diff);
                maxDiff = std::max(maxDiff, diff);
            }
        }
        if (nTrials > 0) {
            meanDiff /= nTrials;
        }
    }
    
    void RbtChromTest::measureCrossoverDiff(RbtChromElement* chrom,
                                            RbtInt nTrials,
                                            RbtDouble& meanDiff,
                                            RbtDouble& minDiff,
                                            RbtDouble& maxDiff) {
        RbtDoubleList refVec;
        chrom->GetVector(refVec);
        meanDiff = 0.0;
        minDiff = 0.0;
        maxDiff = 0.0;
        RbtChromElementPtr chrom2 = chrom->clone();
        RbtChromElementPtr chrom3 = chrom->clone();
        RbtChromElementPtr chrom4 = chrom->clone();
        for (RbtInt iTrial = 0; iTrial < nTrials; iTrial++) {
            chrom->Randomise();
            chrom2->Randomise();
            Rbt::Crossover(chrom, chrom2, chrom3, chrom4);
            RbtInt i(0);
            RbtDouble diff3 = chrom3->CompareVector(refVec,i);
            i = 0;
            RbtDouble diff4 = chrom4->CompareVector(refVec,i);
            meanDiff += diff3;
            meanDiff += diff4;
            if (iTrial == 0) {
                minDiff = std::min(diff3, diff4);
                maxDiff = std::max(diff3, diff4);
            }
            else {
                minDiff = std::min(minDiff, diff3);
                minDiff = std::min(minDiff, diff4);
                maxDiff = std::max(maxDiff, diff3);
                maxDiff = std::max(maxDiff, diff4);
            }
        }
        if (nTrials > 0) {
            meanDiff /= (nTrials*2);
        }
}
