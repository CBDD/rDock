/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtChromOccupancyRefData.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
