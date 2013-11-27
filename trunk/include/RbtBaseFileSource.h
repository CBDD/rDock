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

//Abstract base class for file-based input
                           
#ifndef _RBTBASEFILESOURCE_H_
#define _RBTBASEFILESOURCE_H_

#include <fstream>
using std::ifstream;

#include "RbtConfig.h"

//useful typedefs
typedef RbtString RbtFileRec;
typedef vector<RbtFileRec> RbtFileRecList;
typedef RbtFileRecList::iterator RbtFileRecListIter;

//Max line length expected in file
const int MAXLINELENGTH = 255;

class RbtBaseFileSource
{
 public:
  //Constructors
  //RbtBaseFileSource(const char* fileName);
  RbtBaseFileSource(const RbtString& fileName);
  RbtBaseFileSource(const RbtString& fileName, const RbtString& strRecDelim);

  //Default destructor
  virtual ~RbtBaseFileSource();

  //Public methods

  RbtString GetFileName();
  //void SetFileName(const char* fileName);
  void SetFileName(const RbtString& fileName);

  //Status and StatusOK parse the file to check for errors
  RbtBool StatusOK();
  RbtError Status();

  //FileStatus and FileStatusOK just try and read the file
  RbtBool FileStatusOK();
  RbtError FileStatus();

  //Multi-record methods
  RbtBool isMultiRecordSupported() {return m_bMultiRec;}
  void NextRecord();
  void Rewind();


 protected:
  //////////////////////////////////////////////////////
  //Protected member functions
  //PURE VIRTUAL - MUST BE OVERRIDDEN IN DERIVED CLASSES
  virtual void Parse() throw (RbtError)= 0;
  // chance to give false when record delimiter at the beginning
  void Read(RbtBool aDelimiterAtEnd=true) throw (RbtError);
  //////////////////////////////////////////////////////

  //Protected data
 protected:
  RbtBool m_bParsedOK;//For use by Parse
  RbtFileRecList m_lineRecs;//Allow direct access to cache from derived classes


  //Private member functions
 private:
  RbtBaseFileSource();//Disable default constructor
  RbtBaseFileSource(const RbtBaseFileSource&);//Copy constructor disabled by default
  RbtBaseFileSource& operator=(const RbtBaseFileSource&);//Copy assignment disabled by default
  void Open() throw (RbtError);
  void Close();
  void ClearCache();

  //Private data
 private:
  RbtString m_strFileName;
  RbtBool m_bReadOK;//For use by Read
  ifstream m_fileIn;
  char* m_szBuf;//Line buffer
  RbtBool m_bFileOpen;//Keep track of whether we've opened the file or not
  RbtBool m_bMultiRec;//Is file multi-record ?
  RbtString m_strRecDelim;//Record delimiter
};

#endif //_RBTBASEFILESOURCE_H_
