/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtBaseBiMolTransform.h#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Base class for bimolecular transforms. Overrides Update()

#ifndef _RBTBASEBIMOLTRANSFORM_H_
#define _RBTBASEBIMOLTRANSFORM_H_

#include "RbtBaseTransform.h"
#include "RbtModel.h"

class RbtBaseBiMolTransform : public RbtBaseTransform
{
	public:
	//Class type string
	static RbtString _CT;
	////////////////////////////////////////
	//Constructors/destructors
	virtual ~RbtBaseBiMolTransform();

	////////////////////////////////////////
	//Public methods
	////////////////

	RbtModelPtr GetReceptor() const;
	RbtModelPtr GetLigand() const;
	RbtModelList GetSolvent() const;
		
	//Override RbtObserver pure virtual
	//Notify observer that subject has changed
	virtual void Update(RbtSubject* theChangedSubject);
	
	
	protected:
	////////////////////////////////////////
	//Protected methods
	///////////////////
 	RbtBaseBiMolTransform(const RbtString& strClass, const RbtString& strName);

	//PURE VIRTUAL - Derived classes must override
	virtual void SetupReceptor() = 0;//Called by Update when receptor is changed
	virtual void SetupLigand() = 0;//Called by Update when ligand is changed
	virtual void SetupSolvent() {};//Called by Update when ligand is changed
	virtual void SetupTransform() = 0;//Called by Update when either model has changed
	
	private:
	////////////////////////////////////////
	//Private methods
	/////////////////
  
	protected:
	////////////////////////////////////////
	//Protected data
	////////////////

      
	private:
	////////////////////////////////////////
	//Private data
	//////////////
	RbtModelPtr m_spReceptor;
	RbtModelPtr m_spLigand;
	RbtModelList m_solventList;
};

#endif //_RBTBASEBIMOLTRANSFORM_H_
