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

#ifndef XA_HEAD1_H
#define XA_HEAD1_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "IXAHWave.h"
#include "XArray1D.h"

namespace xar
{
//---------------------------------------------------------------------------
//	FORWARD REFERENCES
//
//---------------------------------------------------------------------------
//	CONSTANT DEFINITIONS
//
//---------------------------------------------------------------------------
//	MACRO DEFINITIONS
//
//---------------------------------------------------------------------------
//	ENUMERATED DATA TYPES
//
//---------------------------------------------------------------------------
//	STRUCTURE DEFINITIONS
//
//---------------------------------------------------------------------------
//	IN-LINE FUNCTION DEFINITIONS
//
	//
	// The following functions are not strictly necessary, but they can make 
	// template-level programming with XArrays and Waveheads easier
	//

	//! Returns a pointer to the const IXAHWave1D interface or zero if it is unavailable
	template <class T> inline const IXAHWave1D* GetIXAHWave1D(const XArray1D<T>& rXAr1D)
	{
		return dynamic_cast<const IXAHWave1D*>(rXAr1D.GetHeadPtr());
	}

	//! Returns a pointer to the IXAHWave1D interface or zero if it is unavailable
	template <class T> inline IXAHWave1D* GetIXAHWave1D(XArray1D<T>& rXAr1D)
	{
		return dynamic_cast<IXAHWave1D*>(rXAr1D.GetHeadPtr());
	}

	//! Returns X-step or 1 if an IXAHWave1D interface is unavailable
	template <class T> double GetStep(const XArray1D<T>& rXAr1D)
	{
		const IXAHWave1D* ph1 = GetIXAHWave1D(rXAr1D);
		if (!ph1) return 1;
		else return ph1->GetStep(rXAr1D.size());
	}

	//! Returns lower X-boundary or 0 if an IXAHWave1D interface is unavailable
	template <class T> inline double GetXlo(const XArray1D<T>& rXAr1D)
	{
		const IXAHWave1D* ph1 = GetIXAHWave1D(rXAr1D);
		if (!ph1) return 0;
		else return ph1->GetXlo();
	}

	//! Resizes the 1D array and the head as well if it is present
	template <class T> void ResizeH(XArray1D<T>& rXAr1D, index_t iNewDim, T tVal = T())
	{
		IXAHWave1D* ph1 = GetIXAHWave1D(rXAr1D);
		if (ph1) ph1->Resize(rXAr1D.size(), iNewDim);
		rXAr1D.Resize(iNewDim, tVal);
	}

	//! Extracts a sub-XArray1D into another XArray1D with an appropriate head
	template <class T> void GetSubarrayH(const XArray1D<T>& rXAr1DSrc, index_t iBegin, index_t iEnd, XArray1D<T>& rDestSubXArray)
	{
		rXAr1DSrc.GetSubarray(iBegin, iEnd, rDestSubXArray);
		const IXAHWave1D* ph1 = GetIXAHWave1D(rXAr1DSrc);
		if (ph1) 
		{
			IXAHWave1D* ph1new = CreateWavehead1D();
			double xlo = ph1->GetXlo() + GetStep(rXAr1DSrc) * iBegin;
			double xhi = xlo + GetStep(rXAr1DSrc) * (iEnd - iBegin - 1);
			ph1new->SetData(ph1->GetWl(), xlo, xhi);
			rDestSubXArray.SetHeadPtr(ph1new);
		}
	}

	
} // namespace xar closed
//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
#endif	// XA_HEAD1_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
