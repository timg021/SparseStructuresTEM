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

// PMT.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define SUBTRACTION_BASED_METHOD 1

#ifdef SUBTRACTION_BASED_METHOD

#include <chrono>
#include <omp.h>
#include "IXAHWave.h"
#include "XArray2D.h"
#include "XArray3D.h"
#include "XA_data.h"
#include "XA_move3.h"
#include "XA_file.h"



//!!! NOTE that fftw3 and XArray3D have the same notation for the order of the array dimensions. Both use C-style array structure, i.e. the last index changes the fastest, 
//i.e. in XArray3D(dim1, dim2, dim3) the fastest changing dimension is dim3, and in fftw_r2c_3d(n0, n1, n2) the fastest dimension is n2.
//Note also that these dimensions are associated with the following physical coordinates by default: dim1(n0) <-> nz, dim2(n1) <-> ny, dim3(n2) <-> nx.



using namespace xar;


int main()
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

	try
	{
		printf("\nStarting PMT1 program ...");
		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024];
		FILE* ff0 = fopen("PMT.txt", "rt");
		if (!ff0) throw std::exception("Error: cannot open parameter file PMT.txt.");
		fgets(cline, 1024, ff0); // 1st line - comment
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2nd line: Defocus_distance_MIN,MAX,STEP_in_Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading defocus parameters from input parameter file.");
		double zmin = atof(cparam); // minimum defocus in Angstroms 
		double zmax = atof(cparam1); // maximum defocus in Angstroms 
		double zstep = atof(cparam2); // defocus step in Angstroms
		printf("\nDefocus planes parameters: min = %g, max = %g, step = %g", zmin, zmax, zstep);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // background intensity value
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading background intensity parameter from input parameter file.");
		float finten0 = (float)atof(cparam); // minimum defocus in Angstroms 
		printf("\nBackground intensity value = %g", finten0);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // file name base for defocus series of the whole sample
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading defocus series file name base for the whole sample from input parameter file.");
		string filenamebaseIn1 = cparam;
		printf("\nDefocus series file name base for the whole sample = %s", filenamebaseIn1.c_str());
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of different atom types
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of atom types from input parameter file.");
		index_t natomtypes = (index_t)atoi(cparam); 
		printf("\nNumber of different atom types = %zd", natomtypes);
		vector<index_t> natom(natomtypes);
		vector<string> filenamebaseIn2(natomtypes); // file name bases for defocus series of different single atoms
		vector< vector<double> > rpos2(natomtypes); // vector of XYZ positions of template atoms
		for (index_t nat = 0; nat < natomtypes; nat++)
		{
			rpos2[nat].resize(3);
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of atoms of this type and file base name for defocus series of a single atom of this type
			if (sscanf(cline, "%s %s %s %s %s %s", ctitle, cparam, cparam1, cparam2, cparam3, cparam4) != 6) throw std::exception("Error reading atom type %d parameters from input parameter file.", int(nat));
			natom[nat] = index_t(atoi(cparam)); // number of atoms of this type
			rpos2[nat][0] = atof(cparam1); // X-coordinate of template atom no. 'nat'
			rpos2[nat][1] = atof(cparam2); // Y-coordinate of template atom no. 'nat'
			rpos2[nat][2] = atof(cparam3); // Z-coordinate of template atom no. 'nat'
			filenamebaseIn2[nat] = cparam4; // file name base for defocus series of single atom of this type
			printf("\nDefocus series file name base for the single atom no.%zd = %s", nat, filenamebaseIn2[nat].c_str());
		}
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // wavelength in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading wavelength parameter from input parameter file.");
		double wl = atof(cparam); // wavelength in input file units (usually, Angstroms). Unfortunately, it is not saved in the GRD files
		printf("\nWavelength = %g (A)", wl);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // average atom size for masking out in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading atom size parameter from input parameter file.");
		double atomsize = atof(cparam); // atom diameter in physical units
		printf("\nAtom diameter = %g (A)", atomsize);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // average atom trace Z-length for masking out in Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::exception("Error reading atom trace Z-length parameters from input parameter file.");
		double atomlength = atof(cparam); // atom "trace" length in the defocus direction in Angstroms to mask "in" the template atom
		double atomsizeZ0 = atof(cparam1); // atom "trace" length in the defocus direction in Angstroms to mask "out" when searching for atoms of the same type
		double atomsizeZ1 = atof(cparam2); // atom "trace" length in the defocus direction in Angstroms to mask "out" when searching for atoms of the next type
		printf("\nAtom trace lengths in Angstroms: 'in' = %g, 'out_same' = %g, 'out_previous' = %g", atomlength, atomsizeZ0, atomsizeZ1);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // high-frequency bandpath radius
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading low frequency limit from input parameter file.");
		int iHPathRad = index_t(atoi(cparam)); // high-frequency bandpath radius: all (abs.)frequencies lower than this one will be zet to zero
		printf("\nHigh frequency lower threshold = %d", iHPathRad);
		index_t iHPathRad2 = index_t(iHPathRad * iHPathRad);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // output file in Vesta XYZ format for detected atom locations
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name for detected atom locations from input parameter file.");
		string filenameOutXYZ = cparam;
		printf("\nOutput XYZ file name = %s", filenameOutXYZ.c_str());
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // optional auxillary data output mode 
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading parameter for selecting output mode from input parameter file.");
		int iCorrArrayOut = atoi(cparam); // if this parameter is 1, the 1st masked array is output, 2 -> 2nd masked array output, 3 -> 3D correlation output is created.
		printf("\nOptional file output mode = %d", iCorrArrayOut);
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // optional output file name base
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading output file name base for saving auxilliary data.");
		string filenamebaseOut = cparam;
		printf("\nAuxilliary output file name baze = %s", filenamebaseOut.c_str());
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::exception("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
		printf("\nNumber of parallel threads = %d", nThreads);
		if (nThreads < 1)
			throw std::exception("The number of parallel threads in input parameter file should be >= 1.");

		fclose(ff0); // close input parameter file

		//************************************ end reading input parameters from file

		// calculate some useful parameters
		index_t ndefocus = index_t((zmax - zmin) / zstep + 0.5); // number of defocus planes, it determines the number of input files to read
		if (ndefocus < 1) throw std::exception("Error: the number of defocus planes is less than 1.");
		int nz = (int)ndefocus;
		printf("\nNumber of defocus planes = %d.", nz);
		int ny = 4, nx = 4; // nx and ny may be overwritten below by data read from input files
		index_t nangles = 1; // !!! nangles values other than 1 are currently not supported in the code below
		double xmin = 0.0, ymin = 0.0;  // default values - may be overwritten below by data read from input files
		double xstep = 1.0, ystep = 1.0; // default values - may be overwritten below by data read from input files
		int karad = int(atomsize / zstep / 2.0 + 0.5); // atom radius in the number of physical z-step units
		int jarad = int(atomsize / ystep / 2.0 + 0.5); // atom radius in the number of physical y-step units - may be overwritten below by data read from input files
		int iarad = int(atomsize / xstep / 2.0 + 0.5); // atom radius in the number of physical x-step units - may be overwritten below by data read from input files
		int karadt = int(atomlength / zstep / 2.0 + 0.5); // 1/2 length of the template atom image "trace" in the defocus direction to mask "in", in the number of physical z-step units
		int karad0 = int(atomsizeZ0 / zstep / 2.0 + 0.5); // 1/2 length of an atom image "trace" in the defocus direction to mask "out" when searching for atoms of the same type, in the number of physical z-step units
		int karad1 = int(atomsizeZ1 / zstep / 2.0 + 0.5); // 1/2 length of an atom image "trace" in the defocus direction to mask "out" when searching for atoms of the next type, in the number of physical z-step units
		printf("\nAtom trace lengths in grid steps: 'in' = %d, 'out_same' = %d, 'out_previous' = %d", 2 * karadt + 1, 2* karad0 + 1, 2* karad1 + 1);
		index_t natomtotal = 0; // total number of found atoms

		// allocate storage for detected atom positions
		vector< vector< vector<index_t> > > vvvatompos(natomtypes); // positions of all atoms and the "correlation coefficient" for each found atom
		for (index_t nat = 0; nat < natomtypes; nat++)
		{
			vvvatompos[nat].resize(natom[nat]); // vvvatompos[nat] is a nat-size vector of vatompos
			for (index_t na = 0; na < natom[nat]; na++) vvvatompos[nat][na].resize(4); // each vector vvvatompos[nat][na] stores (k,j,i) indexes of the position of one atom
		}

		// make a vector of atom type names (extracted from input 1-atom defocus series file names and used for output only)
		vector<string> strAtomNames(natomtypes); // array of atom type names
		for (index_t nat = 0; nat < natomtypes; nat++)
		{
			index_t ii0 = filenamebaseIn2[nat].rfind("\\") + 1;
			index_t ii1 = filenamebaseIn2[nat].find('.', ii0);
			strAtomNames[nat] = filenamebaseIn2[nat].substr(ii0, ii1 - ii0);
		}
		
		
		//**************************************************** start searching for atom positions

		int nzbbb, nybbb, nxbbb, nzccc, nyccc, nxccc; // dimensions of the trimmed 3D template arrays and of the 3D absolute difference array
		XArray3D<float> aaa, bbb, ccc; // 3D defocus series array, 3D single atom defocus series array and the 3D absolute different array
		XArray3DMove<float> aaamove(aaa), bbbmove(bbb), cccmove(ccc); // associated XArray classes for applying masks to 3D arrays
		XArray2D<float> inten; // 2D array for reading/writing series of z-sections of the 3D arrays from/to disk

		for (index_t nat = 0; nat < natomtypes; nat++) // the cycle over the atom type
		{
			if (natom[nat] == 0) continue; // the case where there are 0 atoms of this type

			// load 3D defocus series array
			printf("\n\nNow searching for atoms type no. %zd (%s) ...", nat + 1, strAtomNames[nat].c_str());
			vector<string> infiles;
			printf("\nReading sample defocus series files %s ...", filenamebaseIn1.c_str());
			FileNames(nangles, ndefocus, filenamebaseIn1, infiles);
			for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
			{
				for (int kk = 0; kk < nz; kk++)
				{
					XArData::ReadFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), wl);
					if (nn == 0 && kk == 0)
					{
						nx = (int)inten.GetDim2();
						ny = (int)inten.GetDim1();
						xmin = GetXlo(inten);
						ymin = GetYlo(inten);
						xstep = GetXStep(inten);
						ystep = GetYStep(inten);
						jarad = int(atomsize / ystep / 2.0 + 0.5);
						iarad = int(atomsize / xstep / 2.0 + 0.5);
						aaa.Resize(nz, ny, nx, 0.0f);
						IXAHWave3D* ph3new = CreateWavehead3D();
						ph3new->SetData(wl, zmin, zmax, ymin, ymin + ystep * ny, xmin, xmin + xstep * nx);
						aaa.SetHeadPtr(ph3new);
					}
					else
					{
						if (inten.GetDim1() != ny) throw std::runtime_error("different ny dimension in input file");
						if (inten.GetDim2() != nx) throw std::runtime_error("different nx dimension in input file");
					}
					for (int jj = 0; jj < ny; jj++)
						for (int ii = 0; ii < nx; ii++)
							aaa[kk][jj][ii] = abs(inten[jj][ii] - finten0);
				}
			}
			printf("\nSize of input images: (nx,ny,nz) = (%d, %d, %d); minimums = (%g, %g, %g); steps = (%g, %g, %g).", nx, ny, nz, xmin, ymin, zmin, xstep, ystep, zstep);

			// optionally mask with zeros the vicinity of locations of previously found atoms of other types
			// we allow the option for not masking the locations of previously found atoms, with the idea that duplicate atoms can be localised and then selected on the basis of the correlation coefficient
			if (karad1 >= 0)
			{
				for (index_t natprev = 0; natprev < nat; natprev++)
					for (index_t na = 0; na < natom[natprev]; na++)
						aaamove.FillCylinderPeriodic(vvvatompos[natprev][na][0], vvvatompos[natprev][na][1], vvvatompos[natprev][na][2], karad1, jarad, iarad, 0.0f);
			}

			// optional auxilliary data output
			if (iCorrArrayOut == 1 && nat == natomtypes - 1) // output the masked 1st input array
			{
				printf("\nWriting masked 1st input array in output files %s ...", filenamebaseOut.c_str());
				XArray2D<float> inten(ny, nx);
				IXAHWave2D* ph2new = CreateWavehead2D();
				inten.SetHeadPtr(ph2new);
				ph2new->SetData(wl, ymin, ymin + ystep * ny, xmin, xmin + xstep * nx);
				FileNames(nangles, ndefocus, filenamebaseOut, infiles);
				for (int nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
				{
					for (int kk = 0; kk < nz; kk++)
					{
						for (int jj = 0; jj < ny; jj++)
							for (int ii = 0; ii < nx; ii++)
								inten[jj][ii] = aaa[kk][jj][ii];
						XArData::WriteFileGRD(inten, infiles[nn * nz + kk].c_str(), xar::eGRDBIN);
					}
				}
			}

			// load 3D template array
			printf("\nReading single atom defocus series files %s ...", filenamebaseIn2[nat].c_str());
			FileNames(nangles, ndefocus, filenamebaseIn2[nat], infiles);
			bbb.Resize(nz, ny, nx, 0.0f);
			for (index_t nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
			{
				for (int kk = 0; kk < nz; kk++)
				{
					XArData::ReadFileGRD(inten, infiles[nn * ndefocus + kk].c_str(), wl);
					if (inten.GetDim1() != ny) throw std::runtime_error("different ny dimension in input file");
					if (inten.GetDim2() != nx) throw std::runtime_error("different nx dimension in input file");
					for (int jj = 0; jj < ny; jj++)
						for (int ii = 0; ii < nx; ii++)
							bbb[kk][jj][ii] = abs(inten[jj][ii] - finten0);
				}
			}

			// find the centre of gravity of the second 3D array, i.e. the position of the template atom
			double integ = 0.0, xpos = 0.0, ypos = 0.0, zpos = 0.0, dtemp;
			for (int kk = 0; kk < nz; kk++)
				for (int jj = 0; jj < ny; jj++)
					for (int ii = 0; ii < nx; ii++)
					{
						dtemp = abs(bbb[kk][jj][ii]);
						integ += dtemp;
						xpos += dtemp * ii;
						ypos += dtemp * jj;
						zpos += dtemp * kk;
					}
			xpos /= integ; ypos /= integ; zpos /= integ;
			index_t ipos2 = index_t(xpos + 0.5), jpos2 = index_t(ypos + 0.5), kpos2 = index_t(zpos + 0.5);
			double xpos2 = xmin + xpos * xstep, ypos2 = ymin + ypos * ystep, zpos2 = zmin + zpos * zstep;
			printf("\nCentre of mass position of single atom array in pixels = (%zd, %zd, %zd), and in physical units = (%g, %g, %g).", ipos2, jpos2, kpos2, xpos2, ypos2, zpos2);
			xpos2 = rpos2[nat][0]; ypos2 = rpos2[nat][1]; zpos2 = rpos2[nat][2];
			ipos2 = index_t((xpos2 - xmin) / xstep + 0.5); jpos2 = index_t((ypos2 - ymin) / ystep + 0.5); kpos2 = index_t((zpos2 - zmin) / zstep + 0.5);
			printf("\n Position of this single atom in the parameter file in pixels = (%zd, %zd, %zd), and in physical units = (%g, %g, %g).", ipos2, jpos2, kpos2, xpos2, ypos2, zpos2);

			// trim all pixels outside atomsize vicinity of the centre of mass of the template 1-atom pattern
			if (kpos2 < karadt || bbb.GetDim1() < 1 + kpos2 + karadt || jpos2 < jarad || bbb.GetDim2() < 1 + jpos2 + jarad || ipos2 < iarad || bbb.GetDim3() < 1 + ipos2 + iarad)
				throw std::runtime_error("atomic size and position parameters are inconsistent in input template files");
			bbbmove.Trim(kpos2 - karadt, bbb.GetDim1() - 1 - kpos2 - karadt, jpos2 - jarad, bbb.GetDim2() - 1 - jpos2 - jarad, ipos2 - iarad, bbb.GetDim3() - 1 - ipos2 - iarad);
			nxbbb = (int)bbb.GetDim3(), nybbb = (int)bbb.GetDim2(), nzbbb = (int)bbb.GetDim1();
			printf("\nDimensions of the 3D template array in pixels are: nx = %d, ny = %d, nz = %d.", nxbbb, nybbb, nzbbb);

			// optional auxilliary data output
			if (iCorrArrayOut == 2 && nat == natomtypes - 1) // output the trimmed template array
			{
				printf("\nWriting masked 2nd input array in output files %s ...", filenamebaseOut.c_str());
				FileNames(nangles, karadt * 2 + 1, filenamebaseOut, infiles);
				inten.Resize(nybbb, nxbbb);
				IXAHWave2D* ph2new = CreateWavehead2D();
				ph2new->SetData(wl, ymin + ystep * (jpos2 - jarad), ymin + ystep * (jpos2 - jarad + bbb.GetDim2()), xmin + xstep * (ipos2 - iarad), xmin + xstep * (ipos2 - iarad + bbb.GetDim3()));
				inten.SetHeadPtr(ph2new);
				for (int nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
				{
					for (int kk = 0; kk < nzbbb; kk++)
					{
						for (int jj = 0; jj < nybbb; jj++)
							for (int ii = 0; ii < nxbbb; ii++)
								inten[jj][ii] = bbb[kk][jj][ii];
						XArData::WriteFileGRD(inten, infiles[nn * nzbbb + kk].c_str(), xar::eGRDBIN);
					}
				} 
			}

			// subtract two 3D arrays, shifting the second array around
			printf("\nSubtracting the template array from the defocus series 3D array ...");

			if (nat == 0)
			{
				nzccc = nz - karadt * 2, nyccc = ny - jarad * 2, nxccc = nx - iarad * 2;
				ccc.Resize(nzccc, nyccc, nxccc, 0.0f);
				printf("\nDimensions of the 3D difference array in pixels are: nx = %d, ny = %d, nz = %d.", nxccc, nyccc, nzccc);
			}

			omp_set_num_threads(nThreads);
			#pragma omp parallel default(none) shared(aaa, bbb, ccc, nzccc, nyccc, nxccc, nzbbb, nybbb, nxbbb)
			{
				float adif, asum, anow, bnow;
				float* paaa, * pbbb;
				#pragma omp for schedule(dynamic) nowait
				for (int k = 0; k < nzccc; k++)
					for (int j = 0; j < nyccc; j++)
						for (int i = 0; i < nxccc; i++)
						{
							adif = asum = 0.0f;
							for (int k1 = 0, kk1 = k; k1 < nzbbb; k1++, kk1++)
								for (int j1 = 0, jj1 = j; j1 < nybbb; j1++, jj1++)
								{
									paaa = &aaa[kk1][jj1][i]; pbbb = &bbb[k1][j1][0];
									for (int i1 = 0, ii1 = i; i1 < nxbbb; i1++, ii1++)
									{
										//adif += abs(aaa[kk1][jj1][ii1] - bbb[k1][j1][i1]);
										//asum += aaa[kk1][jj1][ii1] + bbb[k1][j1][i1];
										anow = *paaa++;
										bnow = *pbbb++;
										adif += abs(anow - bnow);
										asum += anow + bnow;
									}
								}
								ccc[k][j][i] = adif / asum; 
						}
			}

			// optional auxilliary data output
			if (iCorrArrayOut == 3 && nat == 0) // output the 3D difference distribution array
			{
				printf("\nWriting 3D absolute difference array in output files %s ...", filenamebaseOut.c_str());
				FileNames(nangles, nzccc, filenamebaseOut, infiles);
				inten.Resize(nyccc, nxccc);
				IXAHWave2D* ph2new = CreateWavehead2D();
				ph2new->SetData(wl, ymin + ystep * jarad, ymin + ystep * (jarad + nyccc), xmin + xstep * iarad, xmin + xstep * (iarad + nxccc));
				inten.SetHeadPtr(ph2new);
				for (int nn = 0; nn < nangles; nn++) // nangles = 1 is assumed
				{
					for (int kk = 0; kk < nzccc; kk++)
					{
						for (int jj = 0; jj < nyccc; jj++)
							for (int ii = 0; ii < nxccc; ii++)
								inten[jj][ii] = ccc[kk][jj][ii];
						XArData::WriteFileGRD(inten, infiles[nn * nzccc + kk].c_str(), xar::eGRDBIN);
					}
				}
			}

			// find the minimums
			printf("\nFinding minimums in the absolute difference 3D array ...");
			index_t kmin, jmin, imin;
			float cccMax = (float)ccc.Norm(eNormMax);
			for (index_t na = 0; na < natom[nat]; na++)
			{
				natomtotal++;

				float amin = ccc.Min3D(kmin, jmin, imin);

				vvvatompos[nat][na][0] = kmin + (index_t)karadt; // absolute z position of the located atom
				vvvatompos[nat][na][1] = jmin + (index_t)jarad; // absolute y position of the located atom
				vvvatompos[nat][na][2] = imin + (index_t)iarad; // absolute x position of the located atom
				vvvatompos[nat][na][3] = index_t((amin * 1.0e+8) + 0.5); // "correlation coefficient" of the located atom

				double xminA = xmin + vvvatompos[nat][na][2] * xstep;
				double yminA = ymin + vvvatompos[nat][na][1] * ystep;
				double zminA = zmin + vvvatompos[nat][na][0] * zstep;

				printf("\nAtom type %zd, atom number %zd:", nat + 1, na + 1);
				printf("\nShift (i,j,k) of the 2nd array relative to the 1st one producing the minimum absolute difference = (%zd, %zd, %zd).", vvvatompos[nat][na][2], vvvatompos[nat][na][1], vvvatompos[nat][na][0]);
				printf("\nAbsolute position (x,y,z) of the detected atom in physical units = (%g, %g, %g).", xminA, yminA, zminA);
				printf("\nAbsolute difference = %g, 'correlation' coefficient = %g", amin, 1.0f - amin);

				// fill the atomsize vicinity of the found minimum by cccMax values, in order to make possible the search for the next smallest minimum
				if (na < natom[nat] - 1) 
					cccmove.FillCylinderPeriodic(kmin, jmin, imin, karad0, jarad, iarad, cccMax);
			} // end if cycle searching for atoms of the current type 'nat'

		} // end of cycle over different atom types


		// bubble-sort the found locations of atoms of each type separately, in accordance with the decreasing correlation coefficient (previously, increasing z-coordinate)
		// (this is done just for increased convenience of checking the output data manually later on)
		vector<index_t> vtemp(3);
		for (index_t nat = 0; nat < natomtypes; nat++)
		{
			index_t n = natom[nat];
			if (n == 0) continue;
			for (index_t i = 0; i < n - 1; i++)
				// Last i elements are already in place    
				for (index_t j = 0; j < n - i - 1; j++)
					// if (vvvatompos[nat][j][0] > vvvatompos[nat][j + 1][0]) // sorting by z
					if (vvvatompos[nat][j][3] > vvvatompos[nat][j + 1][3]) // sorting by correlation coefficient
						std::swap<vector<index_t> >(vvvatompos[nat][j], vvvatompos[nat][j + 1]);
		}

		// write the locations of the detected atoms into an XYZ file compatible with Vesta XYZ format
		FILE* ff = fopen(filenameOutXYZ.c_str(), "wt");
		printf("\n\nWriting output file %s in Vesta XYZ format ...\n", filenameOutXYZ.c_str());
		fprintf(ff, "%zd\n", natomtotal); // number of detected atoms
		fprintf(ff, "%s\n", "Atom positions detected by PMT"); // free-form file info line
		for (index_t nat = 0; nat < natomtypes; nat++)
			for (index_t na = 0; na < natom[nat]; na++)
				fprintf(ff, "%s %f %f %f %f\n", strAtomNames[nat].c_str(), xmin + vvvatompos[nat][na][2] * xstep, ymin + vvvatompos[nat][na][1] * ystep, zmin + vvvatompos[nat][na][0] * zstep, 1.0 - vvvatompos[nat][na][3] / 1.0e+8);
		fclose(ff);

	}
	catch (std::exception& E)
	{
		printf("\n\n!!!Exception: %s\n", E.what());
	}

	std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
	printf("\n\nMain program finished. Execution time = %I64d s.", std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());

	printf("\nPress any key to exit..."); getchar();
	return 0;

}


#endif //ifdef SUBTRACTION_BASED_METHOD