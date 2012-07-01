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

#include "RbtSiteMapperFactory.h"

#include "RbtSphereSiteMapper.h"
#include "RbtLigandSiteMapper.h"
#include "RbtFileError.h"


//Parameter name which identifies a scoring function definition 
RbtString RbtSiteMapperFactory::_MAPPER("SITE_MAPPER");

////////////////////////////////////////
//Constructors/destructors
RbtSiteMapperFactory::RbtSiteMapperFactory() {

}

RbtSiteMapperFactory::~RbtSiteMapperFactory() {

}

////////////////////////////////////////
//Public methods
////////////////
//Creates a single site mapper object of named class
RbtSiteMapper* RbtSiteMapperFactory::Create(const RbtString& strMapperClass, const RbtString& strName) throw (RbtError) {
  if (strMapperClass == RbtSphereSiteMapper::_CT) return new RbtSphereSiteMapper(strName);
  if (strMapperClass == RbtLigandSiteMapper::_CT) return new RbtLigandSiteMapper(strName);
  
  throw RbtBadArgument(_WHERE_,"Unknown site mapper "+strMapperClass);
}

//Creates a single site mapper object whose class is defined
//as the value of the SITE_MAPPER parameter in the strName section of the
//parameter file source argument. Also sets the site mapper parameters from the
//remaining parameter values in the current section
//Note: the current section is restored to its original value upon exit
RbtSiteMapper* RbtSiteMapperFactory::CreateFromFile(RbtParameterFileSourcePtr spPrmSource, const RbtString& strName) throw (RbtError) {
  RbtString strOrigSection(spPrmSource->GetSection());
  spPrmSource->SetSection(strName);
  if (spPrmSource->isParameterPresent(_MAPPER)) {
    RbtString strMapperClass(spPrmSource->GetParameterValueAsString(_MAPPER));
    //Create new site mapper according to the string value of _MAPPER parameter
    RbtSiteMapper* pSiteMapper = Create(strMapperClass,strName);
    //Set all the mapper parameters from the rest of the parameters listed
    RbtStringList prmList = spPrmSource->GetParameterList();
    for (RbtStringListConstIter prmIter = prmList.begin(); prmIter != prmList.end(); prmIter++) {
      if ((*prmIter) != _MAPPER) {//Skip _SF parameter
	pSiteMapper->SetParameter(*prmIter,spPrmSource->GetParameterValueAsString(*prmIter));
      }
    }
    spPrmSource->SetSection(strOrigSection);
    return pSiteMapper;
  }
  else {
    spPrmSource->SetSection(strOrigSection);
    throw RbtFileMissingParameter(_WHERE_,"Missing "+ _MAPPER +" parameter in section "+strName);
  }
}

