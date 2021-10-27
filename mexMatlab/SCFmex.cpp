#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <iostream>
using namespace std;

inline double fast_exp(double x) {
	double d;
	*(reinterpret_cast<int*>(&d) + 0) = 0;
	*(reinterpret_cast<int*>(&d) + 1) = static_cast<int>(1512775 * x + 1072632447);
	return d;
}

void getRGBXYSeeds(int STEP, int width, int height, int* seedIndices, int* RealNumSeeds)
{
    const bool hexgrid = false;
	int n;
    int xstrips, ystrips;
    int xerr, yerr;
    double xerrperstrip,yerrperstrip;
    int xoff,yoff;
    int x,y;
    int xe,ye;
    int seedx,seedy;
    int i;

    xstrips = (int)(0.5+(double)(width)/(double)(STEP));
    ystrips = (int)(0.5+(double)(height)/(double)(STEP));
    
    xerr = width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = width - STEP*xstrips;}
    yerr = height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = height- STEP*ystrips;}
    
	xerrperstrip = (double)(xerr)/(double)(xstrips);
	yerrperstrip = (double)(yerr)/(double)(ystrips);
    
	xoff = STEP/2;
	yoff = STEP/2;
    
    n = 0;
	for( y = 0; y < ystrips; y++ )
	{
        ye = (y*yerrperstrip);
		for( x = 0; x < xstrips; x++ )
		{
            xe = (x*xerrperstrip);
            seedx = (x*STEP+xoff+xe);
            if(hexgrid){ seedx = x*STEP+(xoff<<(y&0x1))+xe; if(seedx >= width)seedx = width-1; }//for hex grid sampling
            seedy = (y*STEP+yoff+ye);
            i = seedy*width + seedx;
			seedIndices[n] = i;
			n++;
		}
	}
    *RealNumSeeds = n;
}

void PerformSCF(double* rvec, double* gvec, double* bvec, double* kseedsr, double* kseedsg, double* kseedsb, double* kseedsx, double* kseedsy, int width, int height, int RealNumSeeds,double sigma_r, int STEP,double*& sigmar,double*& sigmag,double*& sigmab,int Numiter,double multi)
{
    int x1, y1, x2, y2;
	double r, g, b;
    
    int itr;
    int n;
    int x,y;
    int i;
    int ind;
    int k;
    int sz = width*height;
	const int numk = RealNumSeeds;
	int offset = multi * STEP;
    
    sigmar      = (double*)mxCalloc(numk,sizeof(double));
    sigmag      = (double*)mxCalloc(numk,sizeof(double));
    sigmab      = (double*)mxCalloc(numk,sizeof(double));
    double* inv         = (double*)mxCalloc(numk,sizeof(double));
    double* sigmax      = (double*)mxCalloc(numk,sizeof(double));
    double* sigmay      = (double*)mxCalloc(numk,sizeof(double));

    double distrgb;
    double distxy;
    double dist;
    double SigmaS = 1.0 / (STEP * STEP);
	double SigmaR = 1.0 / (sigma_r * sigma_r * 255 * 255);
 
	for( itr = 0; itr < Numiter; itr++ )
	{
		for( n = 0; n < numk; n++ )
		{
            x1 = kseedsx[n]-offset; if(x1 < 0) x1 = 0;
            y1 = kseedsy[n]-offset; if(y1 < 0) y1 = 0;
            x2 = kseedsx[n]+offset; if(x2 > width)  x2 = width;
            y2 = kseedsy[n]+offset; if(y2 > height) y2 = height;
            
			for( y = y1; y < y2; y++ )
			{
				for( x = x1; x < x2; x++ )
				{
					i = y*width + x;
                    
					r = rvec[i];
					g = gvec[i];
					b = bvec[i];

                    distrgb =	(r - kseedsr[n])*(r - kseedsr[n]) +
                                (g - kseedsg[n])*(g - kseedsg[n]) +
                                (b - kseedsb[n])*(b - kseedsb[n]);
                    distxy  =   (x - kseedsx[n])*(x - kseedsx[n]) +
                                (y - kseedsy[n])*(y - kseedsy[n]);    
                    dist = fast_exp(- distrgb * SigmaR - distxy * SigmaS) + DBL_MIN;

                    sigmar[n] += dist * r;
					sigmag[n] += dist * g;
					sigmab[n] += dist * b;
					sigmax[n] += dist * x;
					sigmay[n] += dist * y;

                    inv[n] += dist;
				}
			
            }
		}
       
        {for (int k = 0; k < numk; k++)
		{
			if (inv[k] <= 0) inv[k] = 1;
		}}
		{for( int k = 0; k < numk; k++ )
		{
			kseedsr[k] = sigmar[k]/inv[k];
			kseedsg[k] = sigmag[k]/inv[k];
			kseedsb[k] = sigmab[k]/inv[k];
			kseedsx[k] = sigmax[k]/inv[k];
			kseedsy[k] = sigmay[k]/inv[k];
		}}
	}

    inv         = (double*)mxCalloc(sz,sizeof(double));
    sigmar      = (double*)mxCalloc(sz,sizeof(double));
    sigmag      = (double*)mxCalloc(sz,sizeof(double));
    sigmab      = (double*)mxCalloc(sz,sizeof(double));
    for (int n = 0; n < numk; n++)  //1s
	{
         x1 = kseedsx[n]-offset; if(x1 < 0) x1 = 0;
         y1 = kseedsy[n]-offset; if(y1 < 0) y1 = 0;
         x2 = kseedsx[n]+offset; if(x2 > width)  x2 = width;
         y2 = kseedsy[n]+offset; if(y2 > height) y2 = height;

		for (int y = y1; y < y2; y++)
		{
			for (int x = x1; x < x2; x++)
			{
				int i = y*width + x;

				r = rvec[i];
				g = gvec[i];
				b = bvec[i];

                distrgb =	(r - kseedsr[n])*(r - kseedsr[n]) +
                            (g - kseedsg[n])*(g - kseedsg[n]) +
                            (b - kseedsb[n])*(b - kseedsb[n]);
                distxy  =	(x - kseedsx[n])*(x - kseedsx[n]) +
                          	(y - kseedsy[n])*(y - kseedsy[n]);    
                dist = fast_exp(- distrgb * SigmaR - distxy *SigmaS) + DBL_MIN;
                
                sigmar[i] += dist * kseedsr[n];
				sigmag[i] += dist * kseedsg[n];
				sigmab[i] += dist * kseedsb[n];

				inv[i] += dist;
			}
		}
	}

    for (int i = 0; i < sz; i++) {
		sigmar[i] = sigmar[i] / inv[i];
        sigmag[i] = sigmag[i] / inv[i];
        sigmab[i] = sigmab[i] / inv[i];
	}

    mxFree(sigmax);
    mxFree(sigmay);
    mxFree(inv);
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1) {
        mexErrMsgTxt("At least one argument is required.") ;
    } else if(nrhs > 5) {
        mexErrMsgTxt("Too many input arguments.");
    }

    //---------------------------
    // Variable declarations
    //---------------------------
    int SetNumSeeds = 2000;//default value
    double sigma_r = 0.2;
    int Numiter = 4; //The real count is Numiter+1.
    double multi = 2;
    int width;
    int height;
    int sz;
    int i, ii;
    int x, y;
    double* rvec; double* gvec; double* bvec;
    int step;
    int* seedIndices;
    int RealNumSeeds;
    double* kseedsx;double* kseedsy;
    double* kseedsr;double* kseedsg;double* kseedsb;
    int k;
    const mwSize* dims;//int* dims;
    
    unsigned char* imgbytes;
    //---------------------------
    int numelements   = mxGetNumberOfElements(prhs[0]) ;
    mwSize numdims = mxGetNumberOfDimensions(prhs[0]) ;
    dims  = mxGetDimensions(prhs[0]) ;
    imgbytes  = (unsigned char*)mxGetData(prhs[0]) ;//mxGetData returns g void pointer, so cast it
    width = dims[1]; height = dims[0];//Note: first dimension provided is height and second is width
    sz = width*height;
    //---------------------------
    SetNumSeeds  = mxGetScalar(prhs[1]);
    sigma_r = mxGetScalar(prhs[2]);
    Numiter = mxGetScalar(prhs[3])-1;
    multi = mxGetScalar(prhs[4]);
    //---------------------------
    // Allocate memory
    //---------------------------
    rvec    = (double*)mxMalloc( sizeof(double)      * sz ) ;
    gvec    = (double*)mxMalloc( sizeof(double)      * sz ) ;
    bvec    = (double*)mxMalloc( sizeof(double)      * sz ) ;
    seedIndices = (int*)mxMalloc( sizeof(int)     * sz );
    
    if(numelements/sz == 1)//if it is g grayscale image, copy the values directly into the rgb vectors
    {
        for(x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            for(y = 0; y < height; y++)
            {
                i = y*width+x;
                rvec[i] = imgbytes[ii];
                gvec[i] = imgbytes[ii];
                bvec[i] = imgbytes[ii];
                ii++;
            }
        }
    }
    else{
        for(x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
        {
            for(y = 0; y < height; y++)
            {
                i = y*width+x;
                rvec[i] = imgbytes[ii];
                gvec[i] = imgbytes[ii+sz];
                bvec[i] = imgbytes[ii+sz+sz];
                ii++;
            }
        }
    }


    //---------------------------
    // Find seeds
    //---------------------------
    step = sqrt((double)(sz)/(double)(SetNumSeeds));
    getRGBXYSeeds(step,width,height,seedIndices,&RealNumSeeds);
    
    kseedsx    = (double*)mxMalloc( sizeof(double)      * RealNumSeeds ) ;
    kseedsy    = (double*)mxMalloc( sizeof(double)      * RealNumSeeds ) ;
    kseedsr    = (double*)mxMalloc( sizeof(double)      * RealNumSeeds ) ;
    kseedsg    = (double*)mxMalloc( sizeof(double)      * RealNumSeeds ) ;
    kseedsb    = (double*)mxMalloc( sizeof(double)      * RealNumSeeds ) ;
    for(k = 0; k < RealNumSeeds; k++)
    {
        kseedsx[k] = seedIndices[k]%width;
        kseedsy[k] = seedIndices[k]/width;
        kseedsr[k] = rvec[seedIndices[k]];
        kseedsg[k] = gvec[seedIndices[k]];
        kseedsb[k] = bvec[seedIndices[k]];
    }


    double* sigmar;
    double* sigmag;
    double* sigmab;
    double* outR;
    double* outG;
    double* outB;
    PerformSCF(rvec, gvec, bvec, kseedsr,kseedsg,kseedsb,kseedsx,kseedsy,width,height,RealNumSeeds,sigma_r,step,sigmar,sigmag,sigmab,Numiter,multi);
    
    plhs[0] = mxCreateNumericMatrix(height,width,mxDOUBLE_CLASS,mxREAL);
    plhs[1] = mxCreateNumericMatrix(height,width,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateNumericMatrix(height,width,mxDOUBLE_CLASS,mxREAL);
    outB = (double*)mxGetData(plhs[2]);
    outG = (double*)mxGetData(plhs[1]);
    outR = (double*)mxGetData(plhs[0]);
    for(x = 0, ii = 0; x < width; x++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
        for(y = 0; y < height; y++)
        {
            i = y*width+x;
            outR[ii] = sigmar[i];
            outG[ii] = sigmag[i];
            outB[ii] = sigmab[i];
            ii++;
        }
    }

    //---------------------------
    // Deallocate memory
    //---------------------------
    mxFree(rvec);
    mxFree(gvec);
    mxFree(bvec);
    mxFree(seedIndices);
    mxFree(kseedsx);
    mxFree(kseedsy);
    mxFree(kseedsr);
    mxFree(kseedsg);
    mxFree(kseedsb);
    mxFree(sigmar);
    mxFree(sigmag);
    mxFree(sigmab);
}
