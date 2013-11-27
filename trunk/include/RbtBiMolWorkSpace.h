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

//Subclass of RbtWorkSpace for simulations involving
//Receptor = Model #0
//Ligand   = Model #1
//Solvent  = Model #2 upwards
#ifndef _RBTBIMOLWORKSPACE_H_
#define _RBTBIMOLWORKSPACE_H_

#include "RbtWorkSpace.h"
#include "RbtPopulation.h"

class RbtBiMolWorkSpace : public RbtWorkSpace
{
	enum eModelID {
		RECEPTOR = 0,
		LIGAND = 1,
		SOLVENT = 2 //Flexible waters are represented as separate models, from 2
	};
	
	public:
		////////////////////////////////////////
		//Constructors/destructors
		
		RbtBiMolWorkSpace();
		virtual ~RbtBiMolWorkSpace();
	
		////////////////////////////////////////
		//Public methods
		////////////////
		
		//Model handling
		RbtModelPtr GetReceptor() const;
		RbtModelPtr GetLigand() const;
		RbtModelList GetSolvent() const;
		RbtBool hasSolvent() const;
		void SetReceptor(RbtModelPtr spReceptor);
		void SetLigand(RbtModelPtr spLigand);
		void SetSolvent(RbtModelList solventList);
		void RemoveSolvent();
        void UpdateModelCoordsFromChromRecords(RbtBaseMolecularFileSource* pSource, RbtInt iTrace);
        
		//Model I/O
		//Saves ligand to file sink
		virtual void Save(RbtBool bSaveScores=true);
		virtual void SaveHistory(RbtBool bSaveScores=true);

	protected:
		////////////////////////////////////////
		//Protected methods
		///////////////////
		
		
	private:
		////////////////////////////////////////
		//Private methods
		/////////////////
		
		RbtBiMolWorkSpace(const RbtBiMolWorkSpace&);//Copy constructor disabled by default		
		RbtBiMolWorkSpace& operator=(const RbtBiMolWorkSpace&);//Copy assignment disabled by default
		void SaveLigand(RbtMolecularFileSinkPtr spSink, RbtBool bSaveScores=true);
		
	protected:
		////////////////////////////////////////
		//Protected data
		////////////////

		
	private:
		////////////////////////////////////////
		//Private data
		//////////////
};

//Useful typedefs
typedef SmartPtr<RbtBiMolWorkSpace> RbtBiMolWorkSpacePtr;//Smart pointer

#endif //_RBTBIMOLWORKSPACE_H_
