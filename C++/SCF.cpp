#include "stdafx.h"
#include <cfloat>
#include <iostream>
#include <fstream>
#include "SCF.h"

inline double fast_exp(double x) {
	double d;
	*(reinterpret_cast<int*>(&d) + 0) = 0;
	*(reinterpret_cast<int*>(&d) + 1) = static_cast<int>(1512775 * x + 1072632447);
	return d;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SCF::SCF()
{
	m_rvec = NULL;
	m_gvec = NULL;
	m_bvec = NULL;

}

SCF::~SCF()
{
	if (m_rvec) delete[] m_rvec;
	if (m_gvec) delete[] m_gvec;
	if (m_bvec) delete[] m_bvec;

}

//==============================================================================
///	DetectLabEdges
//==============================================================================
void SCF::DetectRgbEdges(
	const double*				rvec,
	const double*				gvec,
	const double*				bvec,
	const int&					width,
	const int&					height,
	double*&					edges)
{
	int sz = width*height;

	edges = (double*)calloc(sz, sizeof(double));
	for (int j = 1; j < height - 1; ++j)
	{
		for (int k = 1; k < width - 1; ++k)
		{
			int i = j*width + k;

			double dx = (rvec[i - 1] - rvec[i + 1])*(rvec[i - 1] - rvec[i + 1]) +
				(gvec[i - 1] - gvec[i + 1])*(gvec[i - 1] - gvec[i + 1]) +
				(bvec[i - 1] - bvec[i + 1])*(bvec[i - 1] - bvec[i + 1]);

			double dy = (rvec[i - width] - rvec[i + width])*(rvec[i - width] - rvec[i + width]) +
				(gvec[i - width] - gvec[i + width])*(gvec[i - width] - gvec[i + width]) +
				(bvec[i - width] - bvec[i + width])*(bvec[i - width] - bvec[i + width]);

			edges[i] = (dx + dy);
		}
	}
}

//===========================================================================
///	PerturbSeeds
//===========================================================================
void SCF::PerturbSeeds(
	double*&					kseedsr,
	double*&					kseedsg,
	double*&					kseedsb,
	double*&					kseedsx,
	double*&					kseedsy,
	double*&					edges,
	const int&					numseeds)
{
	const int dx8[8] = { -1, -1,  0,  1, 1, 1, 0, -1 };
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1 };


	for (int n = 0; n < numseeds; ++n)
	{
		int ox = kseedsx[n];//original x
		int oy = kseedsy[n];//original y
		int oind = oy*m_width + ox;

		int storeind = oind;
		for (int i = 0; i < 8; ++i)
		{
			int nx = ox + dx8[i];//new x
			int ny = oy + dy8[i];//new y

			if (nx >= 0 && nx < m_width && ny >= 0 && ny < m_height)
			{
				int nind = ny*m_width + nx;
				if (edges[nind] < edges[storeind])
				{
					storeind = nind;
				}
			}
		}
		if (storeind != oind)
		{
			kseedsx[n] = storeind % m_width;
			kseedsy[n] = storeind / m_width;
			kseedsr[n] = m_rvec[storeind];
			kseedsg[n] = m_gvec[storeind];
			kseedsb[n] = m_bvec[storeind];
		}
	}
}

//===========================================================================
///	GetLABXYSeeds_ForGivenK
//===========================================================================
void SCF::GetRGBXYSeeds_ForGivenK(
	double*&					kseedsr,
	double*&					kseedsg,
	double*&					kseedsb,
	double*&					kseedsx,
	double*&					kseedsy,
	const int&					K,
	const bool&					perturbseeds,
	double*&					edgemag,
	int*						seedIndices,
	int*						numseeds)
{

	const bool hexgrid = false;
	int sz = m_width * m_height;
	int STEP = sqrt(double(sz) / double(K)) + 0.5;
	int i;
	int xoff = STEP / 2;
	int yoff = STEP / 2;

	int xstrips = (int)(0.5 + (double)(m_width) / (double)(STEP));
	int ystrips = (int)(0.5 + (double)(m_height) / (double)(STEP));

	int xerr = m_width - STEP*xstrips; if (xerr < 0) { xstrips--; xerr = m_width - STEP*xstrips; }
	int yerr = m_height - STEP*ystrips; if (yerr < 0) { ystrips--; yerr = m_height - STEP*ystrips; }

	double xerrperstrip = (double)(xerr) / (double)(xstrips);
	double yerrperstrip = (double)(yerr) / (double)(ystrips);
	int xe, ye, seedx, seedy;
	int n(0);
	for (int y = 0; y < ystrips; ++y)
	{
		ye = (y*yerrperstrip);
		for (int x = 0; x < xstrips; ++x)
		{
			xe = (x*xerrperstrip);
			seedx = (x*STEP + xoff + xe);
			if (hexgrid) { seedx = x*STEP + (xoff << (y & 0x1)) + xe; if (seedx >= m_width)seedx = m_width - 1; }//for hex grid sampling
			seedy = (y*STEP + yoff + ye);
			i = seedy*m_width + seedx;
			seedIndices[n] = i;
			++n;
		}
	}
	*numseeds = n;

	kseedsx = new double[n];
	kseedsy = new double[n];
	kseedsr = new double[n];
	kseedsg = new double[n];
	kseedsb = new double[n];

	for (int k = 0; k < n; ++k)
	{
		kseedsx[k] = seedIndices[k] % m_width;
		kseedsy[k] = seedIndices[k] / m_width;
		kseedsr[k] = m_rvec[seedIndices[k]];
		kseedsg[k] = m_gvec[seedIndices[k]];
		kseedsb[k] = m_bvec[seedIndices[k]];
	}

	if (perturbseeds)
	{
		PerturbSeeds(kseedsr, kseedsg, kseedsb, kseedsx, kseedsy, edgemag, n);
	}
}


//===========================================================================
///	PerformSuperpixelSegmentation_VariableSandM
//===========================================================================
void SCF::PerformSoftClusteringFiltering(
	unsigned int*				ubuff,
	double*&					kseedsr,
	double*&					kseedsg,
	double*&					kseedsb,
	double*&					kseedsx,
	double*&					kseedsy,
	const double&				sigr,
	const int&					STEP,
	const int&					numk,
	const int&					NUMITR)
{
	int sz = m_width*m_height;

	int numitr(0);

	//----------------
	int offset = 2 * STEP;
	//----------------

	double*  inv = (double*)calloc(numk, sizeof(double));
	double*  sigmar = (double*)calloc(numk, sizeof(double));
	double*  sigmag = (double*)calloc(numk, sizeof(double));
	double*  sigmab = (double*)calloc(numk, sizeof(double));
	double*  sigmax = (double*)calloc(numk, sizeof(double));
	double*  sigmay = (double*)calloc(numk, sizeof(double));

	double SigmaS(1.0 / (STEP * STEP));
	double SigmaR(1.0 / (sigr * sigr * 255 * 255));

	while (numitr < NUMITR)
	{
		//------
		++numitr;
		//------


		for (int n = 0; n < numk; ++n)
		{
			int y1 = kseedsy[n] - offset; if (y1 < 0) y1 = 0;
			int y2 = kseedsy[n] + offset; if (y2 > m_height) y2 = m_height;
			int x1 = kseedsx[n] - offset; if (x1 < 0) x1 = 0;
			int x2 = kseedsx[n] + offset; if (x2 > m_width)  x2 = m_width;

			for (int y = y1; y < y2; ++y)
			{
				for (int x = x1; x < x2; ++x)
				{
					int i = y*m_width + x;
					
					double r = m_rvec[i];
					double g = m_gvec[i];
					double b = m_bvec[i];

					double distrgb((r - kseedsr[n])*(r - kseedsr[n]) +
						(g - kseedsg[n])*(g - kseedsg[n]) +
						(b - kseedsb[n])*(b - kseedsb[n]));

					double distxy((x - kseedsx[n])*(x - kseedsx[n]) +
						(y - kseedsy[n])*(y - kseedsy[n]));

					//------------------------------------------------------------------------
					double dist = fast_exp(-distrgb * SigmaR - distxy * SigmaS) + DBL_MIN;
					//------------------------------------------------------------------------


					sigmar[n] += dist * m_rvec[i];
					sigmag[n] += dist * m_gvec[i];
					sigmab[n] += dist * m_bvec[i];
					sigmax[n] += dist * x;
					sigmay[n] += dist * y;

					inv[n] += dist;

				}
			}
		}



		{for (int k = 0; k < numk; ++k)
		{
			if (inv[k] <= 0) inv[k] = 1;
		}}
		{for (int k = 0; k < numk; ++k)
		{
			kseedsr[k] = sigmar[k] / inv[k];
			kseedsg[k] = sigmag[k] / inv[k];
			kseedsb[k] = sigmab[k] / inv[k];
			kseedsx[k] = sigmax[k] / inv[k];
			kseedsy[k] = sigmay[k] / inv[k];
		}}


	}

	inv = (double*)calloc(sz, sizeof(double));
	sigmar = (double*)calloc(sz, sizeof(double));
	sigmag = (double*)calloc(sz, sizeof(double));
	sigmab = (double*)calloc(sz, sizeof(double));

	for (int n = 0; n < numk; ++n)
	{
		int y1 = kseedsy[n] - offset; if (y1 < 0) y1 = 0;
		int y2 = kseedsy[n] + offset; if (y2 > m_height) y2 = m_height;
		int x1 = kseedsx[n] - offset; if (x1 < 0) x1 = 0;
		int x2 = kseedsx[n] + offset; if (x2 > m_width)  x2 = m_width;

		for (int y = y1; y < y2; ++y)
		{
			for (int x = x1; x < x2; ++x)
			{
				int i = y*m_width + x;
				
				double r = m_rvec[i];
				double g = m_gvec[i];
				double b = m_bvec[i];

				double distrgb((r - kseedsr[n])*(r - kseedsr[n]) +
					(g - kseedsg[n])*(g - kseedsg[n]) +
					(b - kseedsb[n])*(b - kseedsb[n]));

				double distxy((x - kseedsx[n])*(x - kseedsx[n]) +
					(y - kseedsy[n])*(y - kseedsy[n]));

				//------------------------------------------------------------------------
				double dist = fast_exp(-distrgb * SigmaR - distxy * SigmaS) + DBL_MIN;
				//------------------------------------------------------------------------


				sigmar[i] += dist * kseedsr[n];
				sigmag[i] += dist * kseedsg[n];
				sigmab[i] += dist * kseedsb[n];

				inv[i] += dist;

			}
		}
	}

	for (int i = 0; i < sz; ++i) {
		ubuff[i] = ((int)(sigmar[i] / inv[i]) << 16) + ((int)(sigmag[i] / inv[i]) << 8) + (int)(sigmab[i] / inv[i]);
	}

}

//===========================================================================
///	PerformSCFO_ForGivenK
//===========================================================================
void SCF::PerformSCFO_ForGivenK(
	unsigned int*				ubuff,
	const int					width,
	const int					height,
	const int&					K,//required number of seeds
	const double&				sigr)//weight given to range distance
{
	//--------------------------------------------------
	m_width = width;
	m_height = height;
	int sz = m_width*m_height;
	double* kseedsx; double* kseedsy;
	double* kseedsr; double* kseedsg; double* kseedsb;


	//--------------------------------------------------
	int* seedIndices = new int[sz];
	//--------------------------------------------------
	m_rvec = new double[sz];
	m_gvec = new double[sz];
	m_bvec = new double[sz];
	for (int i = 0; i < sz; ++i)
	{
		m_rvec[i] = ubuff[i] >> 16 & 0xff;
		m_gvec[i] = ubuff[i] >> 8 & 0xff;
		m_bvec[i] = ubuff[i] & 0xff;
	}


	//--------------------------------------------------
	clock_t startTime, endTime;
	startTime = clock();

	bool perturbseeds(false);
	double* edgemag;
	int numseeds;
	if (perturbseeds) DetectRgbEdges(m_rvec, m_gvec, m_bvec, m_width, m_height, edgemag);
	GetRGBXYSeeds_ForGivenK(kseedsr, kseedsg, kseedsb, kseedsx, kseedsy, K, perturbseeds, edgemag, seedIndices, &numseeds);

	int STEP = sqrt(double(sz) / double(K)) + 0.5;//adding a small value in the even the STEP size is too small.
	PerformSoftClusteringFiltering(ubuff, kseedsr, kseedsg, kseedsb, kseedsx, kseedsy, sigr, STEP, numseeds, 4);

	endTime = clock();
	float dfPassTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	CString temp_value = _T("");
	temp_value.Format(_T("%f%c"), dfPassTime, 's');
	AfxMessageBox(temp_value);

	if (perturbseeds) free(edgemag);
	delete[] kseedsr;
	delete[] kseedsg;
	delete[] kseedsb;
	delete[] kseedsx;
	delete[] kseedsy;
}