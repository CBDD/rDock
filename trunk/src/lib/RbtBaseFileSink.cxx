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

#include "RbtBaseFileSink.h"
#include "RbtFileError.h"

using std::ios;
using std::endl;

////////////////////////////////////////
//Constructors/destructors
//RbtBaseFileSink::RbtBaseFileSink(const char* fileName) :
//  m_strFileName(fileName)
//{
//  _RBTOBJECTCOUNTER_CONSTR_("RbtBaseFileSink");
//}

RbtBaseFileSink::RbtBaseFileSink(const RbtString& fileName) :
  m_strFileName(fileName),m_bAppend(false)
{
  _RBTOBJECTCOUNTER_CONSTR_("RbtBaseFileSink");
}

RbtBaseFileSink::~RbtBaseFileSink()
{
  Write();//Just in case there is anything in the cache
  _RBTOBJECTCOUNTER_DESTR_("RbtBaseFileSink");
}


////////////////////////////////////////
//Public methods
////////////////
void RbtBaseFileSink::SetFileName(const RbtString& fileName)
{
  Write();//Just in case there is anything in the cache
  m_strFileName = fileName;
}


RbtError RbtBaseFileSink::Status()
{
  //For file sinks, all we can is try and open the file for writing and see what we catch
  try {
    Open(true);//Open for append, so as not to overwrite the file
    Close();
    //If we get here then everything is fine
    return RbtError();
  }

  //Got an RbtError
  catch (RbtError& error) {
    Close();
    return error;
  }
}


////////////////////////////////////////
//Protected methods
///////////////////
//DM 11 Feb 1999 - add flag to allow the cache to be written without clearing it
void RbtBaseFileSink::Write(RbtBool bClearCache) throw (RbtError)
{
  //Only write the file if there is anything in the cache
  if (isCacheEmpty())
    return;

  try {
    Open(m_bAppend);//DM 06 Apr 1999 - open for append or overwrite, depending on m_bAppend attribute
    for (RbtStringListConstIter iter = m_lineRecs.begin(); iter != m_lineRecs.end(); iter++) {
		// for some reason the << overload is screwed up in some sstream 
		// implementations so it is worth to pay this "pointless" price in conversion
		string delimited((*iter).c_str());
		m_fileOut << delimited << endl;
		//m_fileOut << *iter << endl;
    }
    Close();
    if (bClearCache)
      ClearCache();//Clear the cache so we don't write the file again
  }

  //Got an RbtError
  catch (RbtError& error) {
    Close();
    if (bClearCache)
      ClearCache();
    throw;
  }
}

//Add a complete line to the cache
void RbtBaseFileSink::AddLine(const RbtString& fileRec)
{
  m_lineRecs.push_back(fileRec);
}

//Replace a complete line in the cache
void RbtBaseFileSink::ReplaceLine(const RbtString& fileRec, RbtUInt nRec)
{
  if (nRec < m_lineRecs.size())
    m_lineRecs[nRec] = fileRec;
}

////////////////////////////////////////
//Private methods
/////////////////
void RbtBaseFileSink::Open(RbtBool bAppend) throw (RbtError)
{
  std::_Ios_Openmode openMode = ios_base::out;
  if (bAppend)
    openMode = openMode | ios_base::app;
  m_fileOut.open(m_strFileName.c_str(), openMode);
  if (!m_fileOut)
    throw RbtFileWriteError(_WHERE_,"Error opening "+m_strFileName);
}

void RbtBaseFileSink::Close()
{
  m_fileOut.close();
}

void RbtBaseFileSink::ClearCache()
{
  m_lineRecs.clear();
}
