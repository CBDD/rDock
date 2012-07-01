/*This file is part of Rdock.

    Rdock is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Rdock is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Rdock.  If not, see <http://www.gnu.org/licenses/>.*/

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
