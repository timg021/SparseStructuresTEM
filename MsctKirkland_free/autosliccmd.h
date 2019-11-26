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

//@@@@@@ start TEG code
#include <string>
#include <vector>

#define TEG_MULTITHREADED 1 //if defined, multithreaded execution is used

// TEG class for keeping a count across multiple threads (e.g. for counting the number of active threads)
class Counter_Obj
{
	static int counter;
	static bool isUpdated;
public:
	Counter_Obj() { counter++; }
	~Counter_Obj() { counter--; }
	int GetCount() { return counter;  }
	void SetUpdated(bool status) { isUpdated = status; }
	bool GetUpdated() { return isUpdated; }
	void SetTerminate() { counter = -100000; }
};

int autosliccmd(std::vector<std::string> params, std::vector<double> defocus, std::vector<std::string> fileout);
//@@@@@ end TEG code

