/*=================================================================
 *
 * FASTINTERP2.C    .MEX file for performing fast interpolation
 *
 * y = fastinterp2(x, nbns);
 * 
 * inputs:  x = vector to interpolate
 *          nbns = number of times to subdivide x
 * 
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *x, *y, *nbns, dt, a, b;
  int xlen, ylen, j, i, n, xwid, k, k1, nhlf, nhlf2;

  if (nrhs != 2) 
    mexErrMsgTxt("fastinterp.c needs 2 input args: x, nbns"); 
  
  /*  Read inputs into appropriate local variables */
  xlen = mxGetM(prhs[0]);
  xwid = mxGetN(prhs[0]);
  x = mxGetPr(prhs[0]);
  nbns = mxGetPr(prhs[1]);

  dt = 1/nbns[0];
  n = (int) nbns[0];
  nhlf = (int) n/2;
  nhlf2 = n-nhlf;
  ylen = xlen*n;
  
  /* Create memory for output variable */
  plhs[0] = mxCreateDoubleMatrix(ylen,xwid,mxREAL);
  y = mxGetPr(plhs[0]);
  

  if (xwid == 1)  /* if x has 1 column */
  {
      for (j=0;j<nhlf2;j++)
          y[j] = ((j+1+nhlf)*dt)*x[0];
      
      for (i=1;i<xlen;i++)
      {
          for (j=0;j<n;j++)
              y[i*n-nhlf+j] = dt*((j+1)*x[i] + (n-j-1)*x[i-1]);
      }
      for (j=0;j<nhlf;j++)
          y[ylen-nhlf+j] = ((n-j-1)*dt)*x[xlen-1];
  }
  else  /* if x has multiple cols */
  {
      for (j=0;j<nhlf2;j++)
      {
          a = dt*(j+1+nhlf);
          for (k=0;k<xwid;k++)
              y[j+k*ylen] = a*x[k*xlen];
      }
      for (i=1;i<xlen;i++)
      {
          for (j=0;j<n;j++)
          {
              k1 = i*n+j-nhlf;
              a = dt*(j+1);
              b = dt*(n-j-1);
              for (k=0;k<xwid;k++)
                  y[k1+k*ylen] = a*x[i+k*xlen] + b*x[i-1+k*xlen];
          }
      }
      for (j=0;j<nhlf;j++)
      {
          a = dt*(n-j-1);
          k1 = ylen-nhlf+j;
          for (k=0;k<xwid;k++)
              y[k1+k*ylen] = a*x[xlen-1+k*xlen];
      }
     
      
  }
}

