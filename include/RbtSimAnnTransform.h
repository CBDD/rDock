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

//Simple, non-adaptive simulated annealing protocol
#ifndef _RBTSIMANNTRANSFORM_H_
#define _RBTSIMANNTRANSFORM_H_

#include "RbtBaseBiMolTransform.h"
#include "RbtRand.h"
#include "RbtChromElement.h"

//Simple class to keep track of Monte Carlo sampling statistics
class RbtMCStats {
	public:
	RbtMCStats();
	void Init(RbtDouble score);
	void InitBlock(RbtDouble score);
	void Accumulate(RbtDouble score, RbtBool bAccepted);
	RbtDouble Mean() const;
	RbtDouble Variance() const;
	RbtDouble AccRate() const;
	RbtDouble _total;
	RbtDouble _total2;
	RbtDouble _blockInitial;
	RbtDouble _blockFinal;
	RbtDouble _blockMin;
	RbtDouble _blockMax;
	RbtDouble _initial;
	RbtDouble _final;
	RbtDouble _min;
	RbtDouble _max;
	RbtInt _steps;
	RbtInt _accepted;
};
typedef SmartPtr<RbtMCStats> RbtMCStatsPtr;//Smart pointer



class RbtSimAnnTransform : public RbtBaseBiMolTransform
{
 public:
	//Static data member for class type
	static RbtString _CT;
	//Parameter names
	static RbtString _START_T;
	static RbtString _FINAL_T;
	static RbtString _BLOCK_LENGTH;
	static RbtString _SCALE_CHROM_LENGTH;
	static RbtString _NUM_BLOCKS;
	static RbtString _STEP_SIZE;
	static RbtString _MIN_ACC_RATE;
	static RbtString _PARTITION_DIST;
	static RbtString _PARTITION_FREQ;
	static RbtString _HISTORY_FREQ;
	

  ////////////////////////////////////////
  //Constructors/destructors
 	RbtSimAnnTransform(const RbtString& strName = "SIMANN");
	virtual ~RbtSimAnnTransform();

  ////////////////////////////////////////
  //Public methods
  ////////////////	
  
 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
	virtual void SetupReceptor();//Called by Update when receptor is changed
	virtual void SetupLigand();//Called by Update when ligand is changed
	virtual void SetupTransform();//Called by Update when either model has changed
    void MC(RbtDouble t, RbtInt blockLen, RbtDouble stepSize);
    virtual void Execute();

 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////
  RbtSimAnnTransform(const RbtSimAnnTransform&);//Copy constructor disabled by default
  RbtSimAnnTransform& operator=(const RbtSimAnnTransform&);//Copy assignment disabled by default

 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////

 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtMCStatsPtr m_spStats;
  RbtRand& m_rand;//keep a reference to the singleton random number generator
  RbtRequestPtr m_spPartReq;//Partitioning request
  RbtChromElementPtr m_chrom;//Current chromosome
  RbtDoubleList m_minVector;//Chromosome vector corresponding to overall minimum score
  RbtDoubleList m_lastGoodVector;//Saved chromosome before each MC mutation (to allow revert)
};

//Useful typedefs
typedef SmartPtr<RbtSimAnnTransform> RbtSimAnnTransformPtr;//Smart pointer

#endif //_RBTSIMANNTRANSFORM_H_
