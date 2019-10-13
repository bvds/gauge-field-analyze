#ifndef __FORTRAN_INCLUDED
#define __FORTRAN_INCLUDED
/*
           FORTRAN to C data type matching
*/
typedef double doublereal;
typedef int integer;
typedef void fsub;  // gfortran convension
typedef fsub (*U_fp)(); /* Unknown procedure type (function name) */
typedef /* Subroutine */ fsub (*S_fp)();
typedef doublereal (*D_fp)();
typedef integer ftnlen;
typedef integer logical;
/* Used "nm minresqlpModule.o" to determine this */
#define MINRESQLP __minresqlpmodule_MOD_minresqlp
#define DGEMV dgemv_
#define DNRM2 dnrm2_
#define DAXPY daxpy_

extern int MINRESQLP(
    const integer *n, S_fp aprod, const doublereal *b,
    const doublereal *shift, S_fp msolve, S_fp userOrtho, const logical *disable,
    const integer *nout, const integer *itnlim, const doublereal *rtol,
    const doublereal *abstol,
    const doublereal *maxxnorm, const doublereal *trancond,
    const doublereal *acondLim,
    doublereal *x, integer *istop, integer *itn, doublereal *rnorm,
    doublereal *arnorm, doublereal *xnorm, doublereal *anorm,
    doublereal *acond);
// BLAS2 routine
extern int DGEMV(const char *, const integer *, const integer *,
                 const doublereal *, const doublereal *, const integer *,
                 const doublereal *, const integer *, const doublereal *,
                 doublereal *, const integer *);
extern doublereal DNRM2(const integer *, const doublereal *, const integer *);
extern int DAXPY(const integer *n,  const doublereal *alpha,
                 const doublereal *x, const integer *incx,
                 const doublereal *y, const integer *incy);
#endif
