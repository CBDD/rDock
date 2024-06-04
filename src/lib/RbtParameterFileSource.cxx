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

#include "RbtParameterFileSource.h"

#include <algorithm>  //For std::find

#include "RbtFileError.h"

// Constructors
RbtParameterFileSource::RbtParameterFileSource(const char* fileName): RbtBaseFileSource(fileName) {
    _RBTOBJECTCOUNTER_CONSTR_("RbtParameterFileSource");
}

RbtParameterFileSource::RbtParameterFileSource(const RbtString& fileName): RbtBaseFileSource(fileName) {
    _RBTOBJECTCOUNTER_CONSTR_("RbtParameterFileSource");
}

// Destructor
RbtParameterFileSource::~RbtParameterFileSource() {
    ClearParamsCache();
    _RBTOBJECTCOUNTER_DESTR_("RbtParameterFileSource");
}

////////////////////////////////////////
// Public methods
RbtString RbtParameterFileSource::GetTitle() {
    Parse();
    return m_strTitle;
}

RbtString RbtParameterFileSource::GetVersion() {
    Parse();
    return m_strVersion;
}

// DM 06 June 2000 - limit #parameters to those in current section
RbtUInt RbtParameterFileSource::GetNumParameters() {
    Parse();
    return GetParameterList().size();
}

// DM 06 June 2000 - limit parameters to those in current section
RbtStringList RbtParameterFileSource::GetParameterList() {
    Parse();
    std::vector<std::string> parameters;
    if (current_section != nullptr) {
        for (auto& it : current_section->params)
            parameters.emplace_back(it.second);
    }
    return parameters;
}

RbtVariant RbtParameterFileSource::getCurrentSectionParameter(const std::string& name) {
    if (current_section == nullptr) {
        throw RbtFileSectionNotSet(_WHERE_, "Section not set whent trying to fetch parameter value");
    } else {
        auto params = current_section->params;
        auto iter = params.find(name);
        if (iter == params.end())
            throw RbtFileMissingParameter(_WHERE_, GetFullParameterName(name) + " parameter not found in " + GetFileName());
        else
            return iter->second;
    }
}

// DM 4 Feb 1999 Get a particular named parameter value as a double
RbtDouble RbtParameterFileSource::GetParameterValue(const RbtString& strParamName) {
    Parse();
    return getCurrentSectionParameter(strParamName);
}

// DM 12 Feb 1999 Get a particular named parameter value as a string
RbtString RbtParameterFileSource::GetParameterValueAsString(const RbtString& strParamName) {
    Parse();
    return getCurrentSectionParameter(strParamName);
}

// DM 11 Feb 1999 Check if parameter is present
RbtBool RbtParameterFileSource::isParameterPresent(const RbtString& strParamName) {
    Parse();
    if (current_section == nullptr) {
        throw RbtFileSectionNotSet(_WHERE_, "Section not set whent trying to fetch parameter value");
    } else {
        auto params = current_section->params;
        return params.find(strParamName) != params.end();
    }
}
        

// DM 11 Feb 1999 - section handling
// Parameters can be grouped into named sections
// such that the same parameter name can appear in multiple sections
// Default namespace is the unnamed section.
// Main use is for simulation protocols which may need to define a variable number
// of phases - e.g. high temperature sampling, cooling phase, minimisation phase
// and need the same parameters to appear in each

// Number of named sections
RbtInt RbtParameterFileSource::GetNumSections() {
    Parse();
    return sections.size();
}

// List of section names
RbtStringList RbtParameterFileSource::GetSectionList() {
    Parse();
    std::vector<std::string> section_names;
    for (auto& section : sections) section_names.emplace_back(section.name);
    return section_names;
}

// Get current section name
// This is essentially just a prefix for each parameter name
// so we don't actually need to parse the file to get/set the section name
RbtString RbtParameterFileSource::GetCurrentSectionName() const {
    return (current_section == nullptr) ? "" : current_section->name;
}

// Set current section name
void RbtParameterFileSource::SetCurrentSection(const RbtString& strSection) {
    if (strSection == "")
        current_section = nullptr;
    else
        current_section = &getSectionByName(strSection);
}

RbtParameterFileSource::Section & RbtParameterFileSource::getSectionByName(const std::string& name) const {
    auto it = section_name_mapping.find(name);
    if (it == section_name_mapping.end())
        throw RbtFileUnknownSection(_WHERE_, name);
    else
        return it->second;
}

// Private methods

// Pure virtual in RbtBaseFileSource - needs to be defined here
void RbtParameterFileSource::Parse() {
    const RbtString strRbtKey = "RBT_PARAMETER_FILE_V1.00";
    const RbtString strTitleKey = "TITLE ";
    const RbtString strVersionKey = "VERSION ";
    const RbtString strSectionKey = "SECTION ";
    const RbtString strEndSectionKey = "END_SECTION";

    // Only parse if we haven't already done so
    if (!m_bParsedOK) {
        ClearParamsCache();  // Clear current cache
        Read();              // Read the file

        try {
            RbtFileRecListIter fileIter = m_lineRecs.begin();
            RbtFileRecListIter fileEnd = m_lineRecs.end();

            //////////////////////////////////////////////////////////
            // 1. Check for RBT string on line 1
            if ((*fileIter).find(strRbtKey) != 0)
                throw RbtFileParseError(_WHERE_, "Missing " + strRbtKey + " string in " + GetFileName());

            //////////////////////////////////////////////////////////
            // 2. Parse the rest of file
            while (++fileIter != fileEnd) {
                // Ignore blank lines and comment lines
                if (((*fileIter).length() == 0) || ((*fileIter).at(0) == '#')) {
#ifdef _DEBUG
                    // cout << "Comment record" << endl;
#endif  //_DEBUG
                    continue;
                }
                // Check for Title record
                else if ((*fileIter).find(strTitleKey) == 0) {
                    m_strTitle = *fileIter;
                    m_strTitle.erase(0, strTitleKey.length());
#ifdef _DEBUG
                    // cout << "Title = " << m_strTitle << endl;
#endif  //_DEBUG
                }
                // Check for Version record
                else if ((*fileIter).find(strVersionKey) == 0) {
                    m_strVersion = *fileIter;
                    m_strVersion.erase(0, strVersionKey.length());
#ifdef _DEBUG
                    // cout << "Version = " << m_strVersion << endl;
#endif  //_DEBUG
                }
                // Check for Section record
                else if ((*fileIter).find(strSectionKey) == 0) {
                    RbtString section_name;
                    RbtString strDummy;
                    istringstream istr((*fileIter).c_str());
                    istr >> strDummy >> section_name;
                    if (section_name.empty())
                        throw RbtFileParseError(_WHERE_, "Missing SECTION name in " + GetFileName());
                    else if (section_name_mapping.find(section_name) != section_name_mapping.end())
                        throw RbtFileParseError(_WHERE_, "Duplicate " + section_name + " SECTION name in " + GetFileName());
                    else {
                        sections.emplace_back(section_name);
                        section_name_mapping[sections.back().name] = sections.back();
                        SetCurrentSection(section_name);
                    }
                }
                // Check for End Section record
                else if ((*fileIter).find(strEndSectionKey) == 0) {
                    // Just have to set an empty section name to return to the
                    // default unnamed section
                    SetCurrentSection();
                }
                // Assume everything else is a Key=Parameter name, Value=parameter value pair
                else {
                    RbtString strParamName;
                    // DM 12 May 1999 - read as string then convert to variant
                    RbtString strParamValue;
                    istringstream istr((*fileIter).c_str());
                    istr >> strParamName >> strParamValue;
                    // Prefix the parameter name with the section name and ::
                    // Hopefully, this will ensure unique parameter names between sections
                    std::string fullyQualifiedParamName = GetFullParameterName(strParamName);
                    current_section->params[strParamName] = RbtVariant(strParamValue);
#ifdef _DEBUG
                    // cout << strParamName<< " = " << dParamValue << endl;
#endif  //_DEBUG
                }
            }
            //////////////////////////////////////////////////////////
            // If we get this far everything is OK
            m_bParsedOK = true;
            SetCurrentSection();  // Reset to the unnamed section, in case the final END_SECTION record is missing
        }

        catch (RbtError& error) {
            ClearParamsCache();
            throw;  // Rethrow the RbtError
        }
    }
}

void RbtParameterFileSource::ClearParamsCache() {
    m_strTitle = "";
    m_strVersion = "";
    current_section = nullptr;
    sections.clear();
    section_name_mapping.clear();
}

// Returns the fully qualified parameter name (<section>::<parameter name>)
// Checks if name already contains a section name, if so just returns the name unchanged
// If not, prefixes the name with the current section name
RbtString RbtParameterFileSource::GetFullParameterName(const RbtString& strParamName) {
    // Already has section name present, so return unchanged
    if (strParamName.find("::") != RbtString::npos) return strParamName;
    // Else adorn with current section name
    else
        return GetCurrentSectionName() + "::" + strParamName;
}
