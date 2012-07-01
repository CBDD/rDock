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

#include "RbtPMF.h"

RbtString Rbt::PMFType2Str(RbtPMFType aType)
{
	switch (aType) {
		case CF: return RbtString("CF"); break;
		case CP: return RbtString("CP"); break;
		case cF: return RbtString("cF"); break;
		case cP: return RbtString("cP"); break;
		case C3: return RbtString("C3"); break;
		case CW: return RbtString("CW"); break;
		case CO: return RbtString("CO"); break;
		case CN: return RbtString("CN"); break;
		case C0: return RbtString("C0"); break;
		case NC: return RbtString("NC"); break;
		case NP: return RbtString("NP"); break;
		case NA: return RbtString("NA"); break;
		case ND: return RbtString("ND"); break;
		case NR: return RbtString("NR"); break;
		case N0: return RbtString("N0"); break;
		case NS: return RbtString("NS"); break;
		case OC: return RbtString("OC"); break;
		case OA: return RbtString("OA"); break;
		case OE: return RbtString("OE"); break;
		case OR: return RbtString("OR"); break;
		case OS: return RbtString("OS"); break;
		case OD: return RbtString("OD"); break;
		case OW: return RbtString("OW"); break;
		case P:  return RbtString("P");  break;
		case SA: return RbtString("SA"); break;
		case SD: return RbtString("SD"); break;
		case HL: return RbtString("HL"); break;
		case HH: return RbtString("HH"); break;
		case Zn: return RbtString("Zn"); break;
		case CL: return RbtString("CL"); break;
		case Mn: return RbtString("Mn"); break;
		case Mg: return RbtString("Mg"); break;
		case F:  return RbtString("F");  break;
		case Fe: return RbtString("Fe"); break;
		case Br: return RbtString("Br"); break;
		case V:  return RbtString("V"); break;
		case PMF_UNDEFINED:  return RbtString("PMF_UNDEFINED"); break;
		default: return RbtString("NO SUCH TYPE"); break;
	}
}


RbtPMFType  Rbt::PMFStr2Type(RbtString aStr)
{
	if( !aStr.compare("CF") )	
		return CF;
	else if( !aStr.compare("CP") )
		return CP; 
	else if( !aStr.compare("cF") )
		return cF;
	else if (!aStr.compare("cP") )
		return cP;
	else if (!aStr.compare("C3") )
		return C3;
	else if (!aStr.compare("CW") )
		return CW;
	else if (!aStr.compare("CO") )
		return CO;
	else if (!aStr.compare("CN") )
		return CN;
	else if (!aStr.compare("C0") )
		return C0;
	else if (!aStr.compare("NC") )
		return NC;
	else if (!aStr.compare("NP") )
		return NP;
	else if (!aStr.compare("NA") )
		return NA;
	else if (!aStr.compare("ND") )
		return ND;
	else if (!aStr.compare("NR") )
		return NR;
	else if (!aStr.compare("N0") )
		return N0;
	else if (!aStr.compare("NS") )
		return NS;
	else if (!aStr.compare("OC") )
		return OC;
	else if (!aStr.compare("OA") )
		return OA;
	else if (!aStr.compare("OE") )
		return OE;
	else if (!aStr.compare("OR") )
		return OR;
	else if (!aStr.compare("OS") )
		return OS;
	else if (!aStr.compare("OD") )
		return OD;
	else if (!aStr.compare("OW") )
		return OW;
	else if (!aStr.compare("P_") || !aStr.compare("P"))
		return P;
	else if (!aStr.compare("SA") )
		return SA;
	else if (!aStr.compare("SD") )
		return SD;
	else if (!aStr.compare("HL") )
		return HL;
	else if (!aStr.compare("HH") )
		return HH;
	else if (!aStr.compare("Zn") )
		return Zn;
	else if (!aStr.compare("CL") )
		return CL;
	else if (!aStr.compare("Mn") )
		return Mn;
	else if (!aStr.compare("Mg") )
		return Mg;
	else if (!aStr.compare("F_") || !aStr.compare("F"))
		return F;
	else if (!aStr.compare("Fe") )
		return Fe;
	else if (!aStr.compare("Br") )
		return Br;
	else if (!aStr.compare("V_") || !aStr.compare("V"))
		return V;
	
	else
		return PMF_UNDEFINED;
}
