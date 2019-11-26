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

#ifndef XA_HEAD2_H
#define XA_HEAD2_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "IXAHWave.h"
#include "XArray2D.h"

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

	//! Returns a pointer to the const IXAHWave2D interface or zero if it is unavailable
	template <class T> inline const IXAHWave2D* GetIXAHWave2D(const XArray2D<T>& rXAr2D)
	{
		return dynamic_cast<const IXAHWave2D*>(rXAr2D.GetHeadPtr());
	}

	//! Returns a pointer to the IXAHWave2D interface or zero if it is unavailable	
	template <class T> inline IXAHWave2D* GetIXAHWave2D(XArray2D<T>& rXAr2D)
	{
		return dynamic_cast<IXAHWave2D*>(rXAr2D.GetHeadPtr());
	}
	
	//! Returns Y-step or 1 if an IXAHWave2D interface is unavailable
	template <class T> inline double GetYStep(const XArray2D<T>& rXAr2D)
	{
		const IXAHWave2D* ph2 = GetIXAHWave2D(rXAr2D);
		if (!ph2) return 1;
		else return ph2->GetYStep(rXAr2D.GetDim1());
	}

	//! Returns X-step or 1 if an IXAHWave2D interface is unavailable
	template <class T> inline double GetXStep(const XArray2D<T>& rXAr2D)
	{
		const IXAHWave2D* ph2 = GetIXAHWave2D(rXAr2D);
		if (!ph2) return 1;
		else return ph2->GetXStep(rXAr2D.GetDim2());
	}

	//! Returns lower Y-boundary or 0 if an IXAHWave2D interface is unavailable
	template <class T> inline double GetYlo(const XArray2D<T>& rXAr2D)
	{
		const IXAHWave2D* ph2 = GetIXAHWave2D(rXAr2D);
		if (!ph2) return 0;
		else return ph2->GetYlo();
	}

	//! Returns lower X-boundary or 0 if an IXAHWave2D interface is unavailable
	template <class T> inline double GetXlo(const XArray2D<T>& rXAr2D)
	{
		const IXAHWave2D* ph2 = GetIXAHWave2D(rXAr2D);
		if (!ph2) return 0;
		else return ph2->GetXlo();
	}

	//! Resizes the 2D array and the head as well if it is present	
	template <class T> inline void ResizeH(XArray2D<T>& rXAr2D, index_t iNewDim1, index_t iNewDim2, T tVal = T())
	{
		IXAHWave2D* ph2 = GetIXAHWave2D(rXAr2D);
		if (ph2) ph2->Resize(rXAr2D.GetDim1(), rXAr2D.GetDim2(), iNewDim1, iNewDim2);
		rXAr2D.Resize(iNewDim1, iNewDim2, tVal);
	}

	//! Extracts a sub-XArray2D into another XArray2D with an appropriate head
	template <class T> void GetSubarrayH(const XArray2D<T>& rXAr2DSrc, index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, XArray2D<T>& rDestSubXArray)
	{
		rXAr2DSrc.GetSubarray(iBeginDim1, iEndDim1, iBeginDim2, iEndDim2, rDestSubXArray);
		const IXAHWave2D* ph2 = GetIXAHWave2D(rXAr2DSrc);
		if (ph2) 
		{
			IXAHWave2D* ph2new = CreateWavehead2D();
			double ylo = ph2->GetYlo() + GetYStep(rXAr2DSrc) * iBeginDim1;
			double yhi = ylo + GetYStep(rXAr2DSrc) * (iEndDim1 - iBeginDim1 - 1);
			double xlo = ph2->GetXlo() + GetXStep(rXAr2DSrc) * iBeginDim2;
			double xhi = xlo + GetXStep(rXAr2DSrc) * (iEndDim2 - iBeginDim2 - 1);
			ph2new->SetData(ph2->GetWl(), ylo, yhi, xlo, xhi);
			rDestSubXArray.SetHeadPtr(ph2new);
		}
	}

	//! Returns a cross-section along X averaged over several rows, with a head when possible
	template <class T> void GetThickXSectionH(const XArray2D<T>& rXAr2DSrc, index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, XArray1D<T>& rDestXArray1D)
	{
		rXAr2DSrc.GetThickXSection(iBeginDim1, iEndDim1, iBeginDim2, iEndDim2, rDestXArray1D);
		const IXAHWave2D* ph2 = GetIXAHWave2D(rXAr2DSrc);
		if (ph2) 
		{
			IXAHWave1D* ph1new = CreateWavehead1D();
			double xlo = ph2->GetXlo() + GetXStep(rXAr2DSrc) * iBeginDim2;
			double xhi = xlo + GetXStep(rXAr2DSrc) * (iEndDim2 - iBeginDim2 - 1);
			ph1new->SetData(ph2->GetWl(), xlo, xhi);
			rDestXArray1D.SetHeadPtr(ph1new);
		}
	}

	//! Returns a cross-section along Y averaged over several columns, with a head when possible
	template <class T> void GetThickYSectionH(const XArray2D<T>& rXAr2DSrc, index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, XArray1D<T>& rDestXArray1D)
	{
		rXAr2DSrc.GetThickYSection(iBeginDim1, iEndDim1, iBeginDim2, iEndDim2, rDestXArray1D);
		const IXAHWave2D* ph2 = GetIXAHWave2D(rXAr2DSrc);
		if (ph2) 
		{
			IXAHWave1D* ph1new = CreateWavehead1D();
			double ylo = ph2->GetYlo() + GetYStep(rXAr2DSrc) * iBeginDim1;
			double yhi = ylo + GetYStep(rXAr2DSrc) * (iEndDim1 - iBeginDim1 - 1);
			ph1new->SetData(ph2->GetWl(), ylo, yhi);
			rDestXArray1D.SetHeadPtr(ph1new);
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
#endif	// XA_HEAD2_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
