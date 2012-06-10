/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtBaseFileSource.h#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
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
