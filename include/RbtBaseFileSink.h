/***********************************************************************
* $Id: //depot/dev/client3/rdock/2006.1/include/RbtBaseFileSink.h#2 $
* Copyright (C) Vernalis (R&D) Ltd 2006
* This file is released under the terms of the End User License Agreement
* in ../docs/EULA.txt
***********************************************************************/

//Abstract base class for text file-based output. Knows how to Open a file,
//Write a text cache to the file, and Close the file.
//It is up to derived classes to provide a Render method to populate the cache
//with text records of the appropriate format.

#ifndef _RBTBASEFILESINK_H_
#define _RBTBASEFILESINK_H_

#include <fstream>
using std::ofstream;

#include "RbtConfig.h"

class RbtBaseFileSink
{
 public:
  ////////////////////////////////////////
  //Constructors/destructors
  //RbtBaseFileSink(const char* fileName);
  RbtBaseFileSink(const RbtString& fileName);

  virtual ~RbtBaseFileSink(); //Default destructor


  ////////////////////////////////////////
  //Public methods
  ////////////////
  RbtString GetFileName() const {return m_strFileName;}
  //void SetFileName(const char* fileName);
  void SetFileName(const RbtString& fileName);
  RbtBool StatusOK() {return Status().isOK();}
  RbtError Status();

  //PURE VIRTUAL - MUST BE OVERRIDDEN IN DERIVED CLASSES
  virtual void Render() throw (RbtError) = 0;

 protected:
  ////////////////////////////////////////
  //Protected methods
  ///////////////////
  //Write the cache to the file
  //DM 11 Feb 1999 - add flag to allow the cache to be written without clearing it
  void Write(RbtBool bClearCache = true) throw (RbtError);
  //Add a complete line to the cache
  void AddLine(const RbtString& fileRec);
  //Replace a complete line in the cache
  void ReplaceLine(const RbtString& fileRec, RbtUInt nRec);
  //Is cache empty
  RbtBool isCacheEmpty() const {return m_lineRecs.empty();}

  //DM 06 Apr 1999 - append attribute is for derived class use only
  RbtBool GetAppend() const {return m_bAppend;}//Get append status (true=append, false=overwrite)
  void SetAppend(RbtBool bAppend) {m_bAppend = bAppend;}//Set append status (true=append, false=overwrite)


 private:
  ////////////////////////////////////////
  //Private methods
  /////////////////

  RbtBaseFileSink(); //Disable default constructor
  RbtBaseFileSink(const RbtBaseFileSink&);//Copy constructor disabled by default
  RbtBaseFileSink& operator=(const RbtBaseFileSink&);//Copy assignment disabled by default

  void Open(RbtBool bAppend=false) throw (RbtError);
  void Close();
  void ClearCache();


 protected:
  ////////////////////////////////////////
  //Protected data
  ////////////////


 private:
  ////////////////////////////////////////
  //Private data
  //////////////
  RbtStringList m_lineRecs;
  RbtString m_strFileName;
  ofstream m_fileOut;
  RbtBool m_bAppend;//If true, Write() appends to file rather than overwriting

};

//Useful typedefs
typedef SmartPtr<RbtBaseFileSink> RbtBaseFileSinkPtr;//Smart pointer

#endif //_RBTBASEFILESINK_H_

