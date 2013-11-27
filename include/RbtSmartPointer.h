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

//Smart pointer template which can handle normal class inheritance
//relationships. Inspired by Omar Loggiodice's ookit SharedPtr and
//David Harvey's "Smart pointer templates in C++" article
//(http://web.ftech.net/~honeyg/articles/smartp.htm)
//
//We don't use a separate reference counting object, but rather store two
//pointers in each smart pointer:
//1) Pointer to underlying object
//2) Pointer to unsigned int (reference counter)
//
//Pros: the only way I could think of to be able to implement assignment of
//      subclass smart pointer to base class smart pointer
//      Possibly faster access to underlying object (only one dereference
//      rather than two)
//Cons: Less memory efficient (store two pointers per smart pointer, not one)
//
//To do: Should probably throw some other exception than std::exception if
//       assert fails
//NOTE: it's a BAD idea to pass a NULL pointer to the SmartPtr<T> constructor

#ifndef _RBTSMARTPOINTER_H_
#define _RBTSMARTPOINTER_H_

//#include "RbtTypes.h"
//#include "RbtError.h"

//Only check smart pointer assertions in debug build
//#ifdef _NDEBUG
//const RbtBool SMART_CHECK = false;
//#else
//const RbtBool SMART_CHECK = true;
//#endif //_NDEBUG


template <class T>
class SmartPtr
{
	public:
	///////////////////////////////////////////////////
	//CONSTRUCTORS, DESTRUCTORS, ASSIGNMENT

	//Default constructor, initialise both pointers to 0
	SmartPtr() : m_pT(NULL),m_pCount(NULL) {};

	//Parameterised constructor
	//Create new counter, initialise to 1
	SmartPtr(T* pT) : m_pT(pT),m_pCount(new RbtUInt(1)){};

	//Copy constructor - copy both pointers, increment counter
	SmartPtr(const SmartPtr<T>& sp) : m_pT(sp.m_pT),m_pCount(sp.m_pCount) {
		if (!Null())
			GetRef();
	};

	//Copy constructor - template version for creating a base class smart pointer from
	//a subclass smart pointer. Need to use public accessor methods for copying the pointers
	template <class T2>
	SmartPtr(const SmartPtr<T2>& sp) {
	  m_pCount = sp.GetCountPtr();
	  T2* pT2 = const_cast<T2*>(sp.Ptr());
	  m_pT = dynamic_cast<T*>(pT2);
	  if (!Null())
	    GetRef();
	}

	//Destructor - decrement counter, delete underlying object
	//if no refs remaining
	~SmartPtr() {
		UnBind();
	};

	//Assignment operator - check for self assignment
	const SmartPtr<T>& operator=(const SmartPtr<T>& sp) {
		if (this != &sp) {
			UnBind();
			m_pCount = sp.m_pCount;
			m_pT = sp.m_pT;
			if (!sp.Null())
				sp.GetRef();
		}
		return *this;
	};

	//Assignment operator - template version for assigning a subclass smart pointer
	//to a base class smart pointer. Will only compile if T2 is a subclass of T
	//No need to worry about self assignment here so can unbind the existing ref
	//first. Need to use public accessor methods to access the pointers on rhs
	template <class T2>
	const SmartPtr<T>& operator=(const SmartPtr<T2>& sp) {
		UnBind();
		m_pCount = sp.GetCountPtr();
		T2* pT2 = const_cast<T2*>(sp.Ptr());
		m_pT = dynamic_cast<T*>(pT2);
		if (!Null())
			GetRef();
		return *this;
	}
	
	/////////////////////////////////////////////////////
	//FRIEND FUNCTIONS - declaration
	//operator== checks for equivalence of underlying pointers
	//friend bool operator==(const SmartPtr<T>& lhs, const SmartPtr<T>& rhs);


	///////////////////////////////////////////////////
	//PUBLIC ACCESSOR METHODS

	//Tests if smart pointer is empty (i.e. does it have a counter?)
	//Doesn't actually tell you if the underlying pointer is null
	bool Null() const {return m_pCount == NULL;};

	//Returns pointer to counter
	unsigned* GetCountPtr() const {return m_pCount;};
	
	//Returns underlying pointer
	T* Ptr() {
		return m_pT;
	};
	const T * Ptr() const {
		return m_pT;
	};

	//Manually remove reference to underlying object
	void SetNull() {
		UnBind();
	};

	//Various means of dereferencing the underlying object
	//Const and non-const versions
	//Throw assert error if smart pointer is empty
	T* operator->() {
		//Assert<RbtAssert>(!SMART_CHECK||!Null());
		return m_pT;
	};
	const T * operator->() const {
		//Assert<RbtAssert>(!SMART_CHECK||!Null());
		return m_pT;
	};
	T& operator*() {
		//Assert<RbtAssert>(!SMART_CHECK||!Null());
		return *m_pT;
	};
	const T& operator*() const {
		//Assert<RbtAssert>(!SMART_CHECK||!Null());
		return *m_pT;
	};
	
	//DM 12 Jun 2000 - yet another way to dereference the smart pointer
	//Allows a smart pointer to be passed as a regular pointer in a function call
	//for example.
	operator T*() const {return m_pT;}
	
	///////////////////////////////////////////////////
	//PRIVATE METHODS AND DATA
	private:
	//Increments counter and returns new value
	unsigned GetRef() const {return ++(*m_pCount);};
	//Decrements counter and returns new value
	//ASSERT: counter should be non-zero before decrementing
	unsigned FreeRef() const {
		//Assert<RbtAssert>(!SMART_CHECK||(*m_pCount)!=0);
		return --(*m_pCount);
	};
	//Decrements counter and deletes underlying object and counter
	//if count is zero
	void UnBind() {
		if (!Null() && FreeRef() == 0)
		{
			delete m_pT;
			delete m_pCount;
		}
		m_pT = NULL;
		m_pCount = NULL;
	};

	T* m_pT;//Pointer to the underlying object
	unsigned* m_pCount;//Pointer to counter
};


/////////////////////////////////////////////////////
//FRIEND FUNCTIONS - implementation
//operator== checks for equivalence of underlying pointers
//DM 22 May 2000 - use .Ptr() method so we don't have to declare as friend
template <class T>
bool operator==(const SmartPtr<T>& lhs, const SmartPtr<T>& rhs)
{
	return lhs.Ptr() == rhs.Ptr();
}

//operator!= checks for non-equivalence of underlying pointers
template <class T>
bool operator!=(const SmartPtr<T>& lhs, const SmartPtr<T>& rhs)
{
	return !(lhs == rhs);
}

//DM 13 Jul 2000 - operator< for sorting containers of smart pointers
template <class T>
bool operator<(const SmartPtr<T>& lhs, const SmartPtr<T>& rhs)
{
	return lhs.Ptr() < rhs.Ptr();
}

//DM 13 Jul 2000 - operator> for sorting containers of smart pointers
template <class T>
bool operator>(const SmartPtr<T>& lhs, const SmartPtr<T>& rhs)
{
	return lhs.Ptr() > rhs.Ptr();
}

#endif //_RBTSMARTPOINTER_H_
