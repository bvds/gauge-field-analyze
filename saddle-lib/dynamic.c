#include <stdlib.h>
#include <stdio.h>
#include "shifts.h"
/*
    Project out infinitesimal gauge transforms from
    the shifts using MINRES/MINRES-QLP
*/

struct {
  SparseRow *matrix;
  unsigned int matrixElements;
  unsigned int matrixDimension;
  doublereal *z;
}  gaugeProductData;

void dynamic(integer n, SparseRow *gauge, integer gaugeDimension,
	     integer gaugeElements, double *in,
	     integer itnlim, doublereal rtol,
	     double *out) {
  int i;
  doublereal *b = malloc(gaugeDimension*sizeof(double));
  doublereal *x = malloc(gaugeDimension*sizeof(double));
  doublereal *z = malloc(n*sizeof(doublereal));
  doublereal shift = 0.0, maxxnorm = 1.0e4;
  doublereal trancond, acondlim = 1.0e15;
  logical disable = 0;
  integer nout = 6, istop, itn;
  doublereal rnorm, arnorm, xnorm, anorm, acond;
  gaugeProductData.matrix = gauge;
  gaugeProductData.matrixElements = gaugeElements;
  gaugeProductData.matrixDimension = n;
  gaugeProductData.z = z;
  
  matrixVector(gaugeDimension, gauge, gaugeElements, in, b);

  trancond = acondlim; // Always use MINRES
  MINRESQLP(
    &gaugeDimension, gaugeProduct, b, &shift, ((void *)0), &disable, &nout,
    &itnlim, &rtol, &maxxnorm, &trancond, &acondlim,
    x, &istop, &itn, &rnorm, &arnorm, &xnorm, &anorm, &acond);

  vectorMatrix(n, gauge, gaugeElements, x, z);

  for(i=0; i<n; i++) {
    out[i] = in[i] - z[i];
  }
  
  free(z);
  free(b);
  free(x);
}

int gaugeProduct(integer *n, doublereal *x, doublereal *y) {
  vectorMatrix(gaugeProductData.matrixDimension, gaugeProductData.matrix,
	       gaugeProductData.matrixElements,
	       x, gaugeProductData.z);
  matrixVector(*n, gaugeProductData.matrix, gaugeProductData.matrixElements,
	       gaugeProductData.z, y);
  return 0;
}

/* in and out must be distinct */
void matrixVector(int n, SparseRow *a, int na, double *in, double *out) {
  int k;
  for(k=0; k<n; k++) {
    out[k] = 0.0;
  }
  for(k=0; k<na; k++) {
    /* printf("Step:  %i %i %i %i %e %e %e\n", n, k, a[k].i, a[k].j, 
       a[k].value, in[a[k].j], out[a[k].i]);
       fflush(stdout); fflush(stderr); */
    out[a[k].i] += a[k].value * in[a[k].j];
  }
}

void vectorMatrix(int n, SparseRow *a, int na, double *in, double *out) {
  int k;
  for(k=0; k<n; k++) {
    out[k] = 0.0;
  }
  for(k=0; k<na; k++) {
    out[a[k].j] += a[k].value * in[a[k].i];
  }
}
