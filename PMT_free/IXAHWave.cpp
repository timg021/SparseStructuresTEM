/*
------------------------------------------------------------------------
Modifications 2019 T.E. Gureyev

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------- NO WARRANTY ------------------
THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
IN NO EVENT SHALL THE AUTHOR BE LIABLE
FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
ANY OTHER PROGRAM). 
------------------------------------------------------------------------
*/

//---------------------------------------------------------------------------
//	INCLUDE FILES
#include "stdafx.h"

#include "XAHWave.h"
//---------------------------------------------------------------------------
//	LOCAL CONSTANT DEFINITIONS
//
//---------------------------------------------------------------------------
//	LOCAL MACRO DEFINITIONS
//
//---------------------------------------------------------------------------
//	LOCAL ENUMERATED DATA TYPES
//
//---------------------------------------------------------------------------
//	LOCAL STRUCTURE DEFINITIONS
//
//---------------------------------------------------------------------------
//	LOCAL IN-LINE FUNCTION DEFINITIONS
//
//---------------------------------------------------------------------------
//	GLOBAL FUNCTION PROTOTYPES
//
//---------------------------------------------------------------------------
//	LOCAL FUNCTION PROTOTYPES
//
//---------------------------------------------------------------------------
//	GLOBAL DATA DEFINITIONS
//
//---------------------------------------------------------------------------
//	STATIC DATA DEFINITIONS
//
/////////////////////////////////////////////////////////////////////////////
//

//! Creates a Wavehead1D object and returns an interface pointer to it
IXAHWave1D* CreateWavehead1D() { return new xar::Wavehead1D(); }

//! Creates a Wavehead2D object and returns an interface pointer to it
IXAHWave2D* CreateWavehead2D() { return new xar::Wavehead2D(); }

//! Creates a Wavehead3D object and returns an interface pointer to it
IXAHWave3D* CreateWavehead3D() { return new xar::Wavehead3D(); }

//! Creates a Wavehead object of the specified type and returns an interface pointer to it
IXAHead* CreateXAHead(const std::string& strHeadType)
{
	std::string temp = strHeadType.substr(0, 10);
	if (!temp.compare("Wavehead1D")) return CreateWavehead1D();
	if (!temp.compare("Wavehead2D")) return CreateWavehead2D();
	if (!temp.compare("Wavehead3D")) return CreateWavehead3D();
	throw std::invalid_argument("invalid_argument 'strHeadType' in CreateXAHead (unknown head type string)"); 
}
/////////////////////////////////////////////////////////////////////////////
//
//	End of Module
//
