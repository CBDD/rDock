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

//Manages the fixed reference data for a model occupancy chromosome element
//A single instance is designed to be shared between all clones of a given element
#ifndef RBTCHROMOCCUPANCYREFDATA_H_
#define RBTCHROMOCCUPANCYREFDATA_H_

#include "RbtConfig.h"

class RbtModel;

class RbtChromOccupancyRefData {
	public:
        //Class type string
        static RbtString _CT;
       //Sole constructor
        RbtChromOccupancyRefData(RbtModel* pModel,
				 RbtDouble stepSize,//mutation step size
                                 RbtDouble threshold);//threshold for enabling/disabling atoms
         virtual ~RbtChromOccupancyRefData();
        
	RbtDouble GetStepSize() const {return m_stepSize;}
        RbtDouble GetModelValue() const;
	void SetModelValue(RbtDouble occupancy);
        RbtDouble GetInitialValue() const {return m_initialValue;}
								
	private:
	RbtModel* m_pModel;
	RbtDouble m_stepSize;
	RbtDouble m_threshold;
	RbtDouble m_initialValue;
};

typedef SmartPtr<RbtChromOccupancyRefData> RbtChromOccupancyRefDataPtr;//Smart pointer

#endif /*RBTCHROMOCCUPANCYREFDATA_H_*/
