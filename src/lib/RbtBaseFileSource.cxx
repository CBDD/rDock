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

#include "RbtBaseFileSource.h"

#include <cstring>

#include "RbtDebug.h"
#include "RbtFileError.h"

// Constructors
// RbtBaseFileSource::RbtBaseFileSource(const char* fileName)
//{
//   m_strFileName = fileName;
//   ClearCache();
//   _RBTOBJECTCOUNTER_CONSTR_("RbtBaseFileSource");
// }

RbtBaseFileSource::RbtBaseFileSource(const RbtString& fileName): m_bFileOpen(false), m_bMultiRec(false) {
    m_strFileName = fileName;
    m_szBuf = new char[MAXLINELENGTH + 1];  // DM 24 Mar - allocate line buffer
    ClearCache();
    _RBTOBJECTCOUNTER_CONSTR_("RbtBaseFileSource");
}

// Multi-record constructor
RbtBaseFileSource::RbtBaseFileSource(const RbtString& fileName, const RbtString& strRecDelim):
    m_bFileOpen(false),
    m_bMultiRec(true),
    m_strRecDelim(strRecDelim) {
    m_strFileName = fileName;
    m_szBuf = new char[MAXLINELENGTH + 1];  // DM 24 Mar - allocate line buffer
    ClearCache();
    _RBTOBJECTCOUNTER_CONSTR_("RbtBaseFileSource");
}

// Default destructor
RbtBaseFileSource::~RbtBaseFileSource() {
    Close();
    ClearCache();
    delete[] m_szBuf;  // DM 24 Mar - delete line buffer
    _RBTOBJECTCOUNTER_DESTR_("RbtBaseFileSource");
}

// Public methods
RbtString RbtBaseFileSource::GetFileName() { return m_strFileName; }

// void RbtBaseFileSource::SetFileName(const char* fileName)
//{
//   Close();
//   ClearCache();
//   m_strFileName = fileName;
// }

void RbtBaseFileSource::SetFileName(const RbtString& fileName) {
    Close();
    ClearCache();
    m_strFileName = fileName;
}

// Status and StatusOK parse the file to check for errors
RbtBool RbtBaseFileSource::StatusOK() { return Status().isOK(); }

RbtError RbtBaseFileSource::Status() {
    // Try parsing the file and see what we catch
    try {
        Parse();
        // If we get here then everything is fine
        return RbtError();
    }

    // Got an RbtError
    catch (RbtError& error) {
        return error;
    }
}

// FileStatus and FileStatusOK just try and read the file
RbtBool RbtBaseFileSource::FileStatusOK() { return FileStatus().isOK(); }

RbtError RbtBaseFileSource::FileStatus() {
    // Try reading the file and see what we catch
    try {
        Read();
        // If we get here then everything is fine
        return RbtError();
    }

    // Got an RbtError
    catch (RbtError& error) {
        return error;
    }
}

// Multi-record methods

// Force the reading of the next record by clearing the cache
// Doesn't actually read the record
void RbtBaseFileSource::NextRecord() {
    if (m_bMultiRec) {
        ClearCache();
    }
}

// Rewind the file back to the first record
void RbtBaseFileSource::Rewind() {
    if (m_bMultiRec) {
        Close();
        ClearCache();
    }
}

// Protected functions

void RbtBaseFileSource::Read(RbtBool delimiterIsAtEnd) {
    if (m_bReadOK) return;  // If we have already read the file, skip
    ClearCache();
    std::string line;
    try {
        Open();
        if (!m_bMultiRec) {
            while (std::getline(m_fileIn, line)) {
                DEBUG_OUT(line << endl);
                m_lineRecs.push_back(line);
            }
            Close();  // original code closed the single record file directly once read
        } else if (delimiterIsAtEnd) {
            while ((std::getline(m_fileIn, line)) && (line != m_strRecDelim)) {
                DEBUG_OUT(line << endl);
                m_lineRecs.push_back(line);
            }
        } else {
            // skip to the header stuff until the first record and the first delimiter line
            while ((std::getline(m_fileIn, line)) && (line != m_strRecDelim)) {
            };
            while ((std::getline(m_fileIn, line)) && (line != m_strRecDelim)) {
                DEBUG_OUT(line << endl);
                m_lineRecs.push_back(line);
            }
        }
        if (m_lineRecs.empty()) throw RbtFileReadError(_WHERE_, "End of file/empty record in " + m_strFileName);
    } catch (RbtError& error) {
        Close();
        throw;
    }
    m_bReadOK = true;
}

// Private functions
void RbtBaseFileSource::Open() {
    // DM 23 Mar 1999 - check if file is already open, to allow Open() to be called redundantly
    if (!m_bFileOpen) m_fileIn.open(m_strFileName.c_str(), ios_base::in);

    // If file did not open, throw an error
    if (!m_fileIn)
        throw RbtFileReadError(_WHERE_, "Error opening " + m_strFileName);
    else
        m_bFileOpen = true;
}

void RbtBaseFileSource::Close() {
    m_fileIn.close();
    m_bFileOpen = false;
}

void RbtBaseFileSource::ClearCache() {
    m_lineRecs.clear();   // Get rid of the previous file records
    m_bReadOK = false;    // Indicate the cache is invalid
    m_bParsedOK = false;  // Tell the Parse function in derived classes that
    // it will have to reparse the file
}
