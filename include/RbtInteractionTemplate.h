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

#ifndef _RBTINTERACTIONTEMPLATE_H_
#define _RBTINTERACTIONTEMPLATE_H_
#include "RbtConfig.h"
#include "RbtBaseGrid.h"

template <class T> class RbtInteractionNode {
	T	property;
	
	public:
		RbtInteractionNode() { _RBTOBJECTCOUNTER_CONSTR_("RbtInteractionNode"); };
		RbtInteractionNode(T aProperty) {
			_RBTOBJECTCOUNTER_CONSTR_("RbtInteractionNode");
			property = aProperty;
		};
		~RbtInteractionNode() { _RBTOBJECTCOUNTER_DESTR_("RbtInteractionNode"); };
		
		T		GetProperty() const { return property; };
		void	SetProperty(T aProperty) { property = aProperty; };
};

template<typename T>	// instead RbtInteractionCenterList
class RbtInteractionNodeList: public std::vector< RbtInteractionNode<T>* >
{ 
}; // std:vector ctors

template<typename T>	// instead RbtInteractionListMap
class RbtInteractionNodeListMap: public std::vector<RbtInteractionNodeList<T> >
{ }; // std:vector ctors 

template <typename T> 
class RbtInteractionNodeListMapIter: public RbtInteractionNodeListMap<T>::iterator 
{ };

template <typename T> 
class RbtInteractionNodeListMapConstIter: public RbtInteractionNodeListMap<T>::const_iterator 
{ };
// Most of the stuff is copied shamelessly from RbtInteractionGrid.[cxx,h]
template <class T> class RbtInteractionGridTemplate : public RbtBaseGrid {
	public:
		static RbtString _CT;
		
  		RbtInteractionGridTemplate(const RbtCoord& gridMin, const RbtCoord& gridStep,
				RbtUInt NX, RbtUInt NY, RbtUInt NZ, RbtUInt NPad=0) {
			_CT = "RbtInteractionGridTemplate";
			CreateMap();
			_RBTOBJECTCOUNTER_CONSTR_("RbtInteractionGridTemplate");
		};
		
		//Constructor reading params from binary stream
		RbtInteractionGridTemplate(istream& istr) {
			_CT = "RbtInteractionGridTemplate";
			CreateMap();
			OwnRead(istr);
			_RBTOBJECTCOUNTER_CONSTR_("RbtInteractionGridTemplate");
		};
		
		~RbtInteractionGridTemplate() {
			ClearInteractionLists();
			_RBTOBJECTCOUNTER_DESTR_("RbtInteractionGridTemplate");
		};
		
		//Copy constructor
		RbtInteractionGridTemplate(const RbtInteractionGridTemplate& aGrid) : RbtBaseGrid(aGrid) {
			CreateMap();
			CopyGrid(aGrid);
			_RBTOBJECTCOUNTER_COPYCONSTR_("RbtInteractionGridTemplate");
		};

		//Copy constructor taking base class argument
		RbtInteractionGridTemplate(const RbtBaseGrid&) {
			CreateMap();
			_RBTOBJECTCOUNTER_COPYCONSTR_("RbtInteractionGridTemplate");
		};

		//Copy assignment
		RbtInteractionGridTemplate& operator=(const RbtInteractionGridTemplate& aGrid) {
			if( this != &aGrid ) {
				ClearInteractionLists();
				RbtBaseGrid::Operator=(aGrid);
				CopyGrid(aGrid);
			}
			return *this;
		};

		//Copy assignment taking base class argument
		RbtInteractionGridTemplate& operator=(const RbtBaseGrid& aGrid) {
			if( this != &aGrid) {
				ClearInteractionLists();
				RbtBaseGrid::Operator=(aGrid);
			}
			return *this;
		};

		virtual void Print(ostream& ostr) const {//Text output
			RbtBaseGrid::Print(ostr);
			OwnPrint(ostr);
		};
		virtual void Write(ostream& ostr) const { //Binary output (serialisation)
			RbtBaseGrid::Write(ostr);
			OwnWrite(ostr);
		};
		virtual void Read(istream& istr) { //Binary input, replaces existing grid
			ClearInteractionLists();
			RbtBaseGrid::Read(istr);
			OwnRead(istr);
		}

		const RbtInteractionNodeList<T>& GetInteractionList(RbtUInt iXYZ) const {
			if (isValid(iXYZ)) {
				return m_intnMap[iXYZ];
			} else {
				return m_emptyList;
			}
		};
		
		const RbtInteractionNodeList<T>& GetInteractionList(const RbtCoord& c) const {
			if (isValid(c)) {
				return m_intnMap[GetIXYZ(c)];
			} else {
				//cout << _CT << "::GetInteractionList," << c << " is off grid" << endl;
				return m_emptyList;
			}
		};

  		void SetInteractionLists(T* pIntn, RbtDouble radius) {
			RbtAtom*	pAtom = pIntn->GetProperty().atom;
			if(NULL == pAtom)
				return;
			const RbtCoord& c = pAtom->GetCoords();
			RbtUIntList sphereIndices;
			GetSphereIndices(c,radius,sphereIndices);
			for (RbtUIntListConstIter sphereIter = sphereIndices.begin(); 
					sphereIter != sphereIndices.end(); 
					sphereIter++) {
				m_intnMap[*sphereIter].push_back(pIntn);
			}
		};

		void ClearInteractionLists() {
			for (RbtInteractionNodeListMap<T>::iterator iter = m_intnMap.begin();  
					iter != m_intnMap.end(); 
					iter++) {
				(*iter).clear();
			}
		};
	
	protected:
		void OwnPrint(ostream& ostr) const {
			ostr << endl << "Class\t" << _CT << endl;
			ostr << "No. of entries in the map: " << m_intnMap.size() << endl;
		};
		void OwnWrite(ostream& ostr) const {};
		void OwnRead(istream& istr) throw (RbtError) {
			CreateMap();
		};

	private:
		RbtInteractionGridTemplate();

		void CopyGrid(const RbtInteractionGridTemplate<T>& aGrid) {
			m_intnMap = aGrid.m_intnMap;
		};
		
		void CreateMap() {
			m_intnMap = RbtInteractionNodeListMap<T>(GetN());
		};

		RbtInteractionNodeListMap<T> m_intnMap;//Used to store the interaction center lists at each grid point
		const RbtInteractionNodeList<T> m_emptyList;//Dummy list used by GetAtomList
};

template <typename T>	// instead RbtInteractionGridPtr
class RbtInteractionGridTemplatePtr : public SmartPtr< RbtInteractionGridTemplate<T> >
{};
#endif //  _RBTINTERACTIONTEMPLATE_H_
