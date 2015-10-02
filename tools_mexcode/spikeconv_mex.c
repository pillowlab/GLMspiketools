/*=================================================================
 *
 * SPIKECONV_MEX.C    .MEX file for performing convolution of the h current 
 *   with a vector of spike times. 
 *   (adds h current beginning in the bin after each spike).  
 *
 * spikeconv_mex(spikeinds, h_currents, twin);
 *
 * Inputs:  spikeinds = vector of spike arrival times (integers)
 *          h_currents = each column is an h current following a spike
 *          twin = [t0 t1], time window over which to compute net h
 *              t0 = index of first bin to compute, t1 = last bin
 *              output has length t1-t0+1
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *spinds, *h, *twin, *hcurrent;
  int nsp, spwid, hlen, hwid, j, i, k, isp, i1, imx, ycol, hcol, rlen, t0, t1, t2;

  if (nrhs != 3) { 
    mexErrMsgTxt("needs 3 input args:  spikeconv_mex(spinds, h_currents, twin)"); 
 }

  /*  Read inputs into appropriate local variables */
  nsp = mxGetM(prhs[0]);
  spwid = mxGetN(prhs[0]);
     if (spwid > 1){
         mexErrMsgTxt("spikeconv_mex: Spike Times Must Be Passed in As a COLUMN Vector!!");
     }
  spinds = mxGetPr(prhs[0]);
  hlen = mxGetM(prhs[1]);
  hwid = mxGetN(prhs[1]);
  h = mxGetPr(prhs[1]);
  twin = mxGetPr(prhs[2]);
  
  /* Create memory for output variable */
  t0 = (int) twin[0]-hlen+1;
  t1 = (int) twin[0]-1;
  t2 = (int) twin[1];
  rlen = t2-t1;
  plhs[0] = mxCreateDoubleMatrix(rlen,hwid,mxREAL);
  hcurrent = mxGetPr(plhs[0]);

  /* Find the first spike times to start adding */
  j = 0;  
  while ( (j<nsp) && (spinds[j]<t0) )
    j++;

  if (j<nsp)
    {
      /* Loop over relevant spike times */
      for ( ; (j<nsp) && (spinds[j]<t2+1) ; j++)
	{
	  isp = (int) spinds[j]-1; /* index of spike */
	  
	  i1 = isp;  /* starting index */
	  if (i1 < t1)
	    i1 = t1;
	  
	  imx = isp+hlen;  /* stopping index */
	  if (imx > t2)
	    imx = t2;

	  /* loop over columns */
	  for (k=0;k<hwid;k++)
	    {
	      ycol = k*rlen-t1;
	      hcol = k*hlen;
	      for (i=i1;i<imx;i++)
		  hcurrent[ycol+i] += h[hcol+i-isp];
	    }
	}
    }
}
