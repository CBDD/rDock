/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtChromPositionRefData.h#3 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Manages the fixed reference data for a position chromosome element
//Also provides methods to map the genotype (COM and Euler angles) onto the
//phenotype (model coords)
//A single instance is designed to be shared between all clones of a given element
#ifndef RBTCHROMPOSITIONREFDATA_H_
#define RBTCHROMPOSITIONREFDATA_H_

#include "RbtModel.h"
#include "RbtDockingSite.h"
#include "RbtEuler.h"
#include "RbtChromElement.h"

class RbtChromPositionRefData {
	public:
        //Class type string
        static RbtString _CT;
        //Reference Cartesian axes
        static const RbtPrincipalAxes CARTESIAN_AXES;
 		RbtChromPositionRefData(const RbtModel* pModel,
                                const RbtDockingSite* pDockSite,
                                RbtDouble transStepSize,//Angstroms
                                RbtDouble rotStepSize,//radians
                                RbtChromElement::eMode transMode=RbtChromElement::FREE,
                                RbtChromElement::eMode rotMode=RbtChromElement::FREE,
                                RbtDouble maxTrans=0.0,//Angstroms
                                RbtDouble maxRot=0.0);//radians
        virtual ~RbtChromPositionRefData();

		RbtInt GetNumStartCoords() const {
            return m_startCoords.size();
        }
        const RbtCoord& GetStartCoord(RbtInt iCoord) const {
            return m_startCoords[iCoord];
        }
        RbtDouble GetTransStepSize() const {return m_transStepSize;}
        RbtDouble GetRotStepSize() const {return m_rotStepSize;}
        RbtChromElement::eMode GetTransMode() const {return m_transMode;}
        RbtChromElement::eMode GetRotMode() const {return m_rotMode;}
        //Chromosome length, excluding FIXED modes (0, 3 or 6)
        RbtInt GetLength() const {return m_length;}
        //Chromosome length for crossover, excluding FIXED modes (0, 1 or 2)
        RbtInt GetXOverLength() const {return m_xOverLength;}
        RbtBool IsTransFixed() const {return m_transMode == RbtChromElement::FIXED;}
        RbtBool IsRotFixed() const {return m_rotMode == RbtChromElement::FIXED;}
        RbtDouble GetMaxTrans() const {return m_maxTrans;}
        RbtDouble GetMaxRot() const {return m_maxRot;}
        const RbtCoord& GetInitialCOM() const {
            return m_initialCom;
        }
        const RbtEuler& GetInitialOrientation() const {
            return m_initialOrientation;
        }
        const RbtQuat& GetInitialQuat() const {
            return m_initialQuat;
        }
        
        void GetModelValue(RbtCoord& com, RbtEuler& orientation) const;
		void SetModelValue(const RbtCoord& com, const RbtEuler& orientation);
								
	private:
        RbtAtomList m_refAtoms;
        RbtAtomRList m_movableAtoms;
        RbtCoordList m_startCoords;
		RbtDouble m_transStepSize;
        RbtDouble m_rotStepSize;
        RbtCoord m_initialCom;
        RbtEuler m_initialOrientation;
        RbtQuat m_initialQuat;
        RbtChromElement::eMode m_transMode;
        RbtChromElement::eMode m_rotMode;
        RbtInt m_length; 
        RbtInt m_xOverLength;
        //Max distance allowed from starting coord
        //Only used if m_transMode == TETHERED
        RbtDouble m_maxTrans;
        //Max rot allowed from starting orientation
         //Only used if m_rotMode == TETHERED
        RbtDouble m_maxRot;
};

typedef SmartPtr<RbtChromPositionRefData> RbtChromPositionRefDataPtr;//Smart pointer

#endif /*RBTCHROMPOSITIONREFDATA_H_*/
