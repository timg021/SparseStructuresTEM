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

#if !defined XARRAYB_H
#define XARRAYB_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_ini.h"

namespace xar
{
//---------------------------------------------------------------------------
//	USING DECLARATIONS
//
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
//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//Class XArrayBase<T>
//
//	Foundation base class for XArray classes
//
/*!
	\brief		Foundation base class template for XArray classes
	\par		Description:
				This is a 'working horse' base class for XArray classes, which handles memory allocation
				and std::vector-type operations. Currently, XArrayBase<T> is just a very slightly augmented
				std::vector<T> class.
	\remarks	Most functionality is currently inherited from std::vector<T>. Still, we need to exlicitly 
				'retranslate' the members of std::vector<T> that are not inherited, but used in the derived classes, 
				i.e. (1) constructors; (2) operator=
	\remarks	Because of the very close link between this class and std::vector, all member functions in this
				class have names starting from small letters (not capital letters, as in all other XArray classes)
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
*/

	template <class T> class XArrayBase : public vector<T>
	{
	// Enumerators
	// Structures
	// Constructors
		//! Constructors are protected to prevent instantiation of objects of this class
	protected:
		//! Default constructor
		XArrayBase() {}
		//! Constructor with a predefined size
		explicit XArrayBase(index_t NumPoints, T tVal = T()) : vector<T>(NumPoints, tVal) {}
		 //! Promotion from a vector
		explicit XArrayBase(const vector<T>& rVector) : vector<T>(rVector) {}
		//! Construction from a raw memory buffer
		XArrayBase(T* ptBufBegin, T* ptBufEnd) : vector<T>(ptBufBegin, ptBufEnd) {} 
		 //! Copy constructor
		XArrayBase(const XArrayBase<T>& rXArrayBase) : vector<T>(rXArrayBase) {}
		//! Destructor
		~XArrayBase() {} 

	// Operators
	public:
		//! Makes this  a (deep) copy of the rXArrayBase
		void operator=(const XArrayBase<T>& rXArrayBase) { vector<T>::operator=(rXArrayBase); }

	// Operations
	public:
		//! Resizes the array to zero and FREES THE MEMORY
		void truncate();
		//! Accepts an external memory buffer with its contents
		// Dirty implementation (depends on the Microsoft implementation of std::vector)
		void acceptMemBuffer(T* ptBufBegin, index_t BufSize); 
		//! Relinquishes the responsibility for the memory area occupied by the internal array
		// Dirty implementation (depends on the Microsoft implementation of std::vector)
		void releaseMemBuffer();

	// Overridables
	public:
	
	// Implementation
	protected:

	private:
	// Member variables	
	// Member functions

	};

}

//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
namespace xar
{
	//! Resizes the array to zero and FREES THE MEMORY
	// after B.Stroustrup, 3rd edition
	template <class T> void XArrayBase<T>::truncate()
	{
		XArrayBase<T> temp(index_t(0)); 
		swap(temp);
	}

	//! Accepts an external memory buffer with its contents
	// This buffer can be legitimately passed to any vector<T>-derived class 
	// without copying(!) by calling the vector::swap function
	// This is a dirty solution for the 'vector accepting an external memory buffer' problem.
	// The solution uses the Visual C++ 6.0 implementaton details of the STL <vector>, namely
	// the protected members _First, _Last and _End. Note also that some related standard
	// features (e.g. 'hint' in allocator.allocate) are in fact NOT IMPLEMENTED in the VC++ 6.0 STL.
/*
	template <class T> void XArrayBase<T>::acceptMemBuffer(T* ptBufBegin, index_t BufSize)
	{
		truncate();
#if _MSC_VER <= 1200 // VC++ 6.0 and lower
	_First = ptBufBegin; _Last = _End = ptBufBegin + BufSize;
#else
	_Myfirst = ptBufBegin; _Mylast = _Myend = ptBufBegin + BufSize;
#endif
	}

	//! Relinquishes the responsibility for the memory area occupied by the internal array	
	// This is a dirty solution for the 'vector surrendering its memory buffer' problem.
	// The solution uses the Visual C++ 6.0 implementaton details of the STL <vector>, namely
	// the protected members _First, _Last and _End. Note also that some related standard
	// features (e.g. 'hint' in allocator.allocate) are in fact NOT IMPLEMENTED in the VC++ 6.0 STL.
	template <class T> void XArrayBase<T>::releaseMemBuffer()
	{
#if _MSC_VER <= 1200 // VC++ 6.0 and lower
		_First = _Last = _End = 0;
#else
		_Myfirst = _Mylast = _Myend = 0;
#endif
	}
*/
}
//---------------------------------------------------------------------------
//	RELATED IN-LINE NON-MEMBER DEFINITIONS
//
//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//
// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArrayBase<char>;
	template class xar::XArrayBase<short>;
	template class xar::XArrayBase<long>;
	template class xar::XArrayBase<float>;
	template class xar::XArrayBase<double>;
	template class xar::XArrayBase<xar::fcomplex>;
	template class xar::XArrayBase<xar::dcomplex>;
#endif

//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
#endif	// XARRAYB_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
