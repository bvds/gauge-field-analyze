/* a C version of simple77.f */
#include <stdio.h>

//MMA.cc
#include "WolframLibrary.h"
#include "WolframCompileLibrary.h"

DLLEXPORT mint WolframLibrary_getVersion(){
  return WolframLibraryVersion;
}
DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData){
  return 0;
}

/* This set of parameters tell TRLAN to compute 5 (NED) smallest
   (LOHI=-1) of a 100 x 100 (NROW) matrix by using a Krylov subspace
   basis size of 30 (MAXLAN).
   MEV defines the size of the arrays used to store computed solutions,
   eval, evec.
*/
#define NROW	100
#define LOHI	-1
#define NED	5
#define MAXLAN	30
#define MEV	10

/* global variables */
int my_pe = 0;

#ifdef CRAY
#define F_TRLAN TRLAN77
#else
#define F_TRLAN trlan77_
#endif

/* prototype of operator (matrix-vector multiplication routine */
void diag_op(int*, int*, double*, int*, double*, int*);
/* prototype of the TRLan77 function (written in FORTRAN)
   NOTE: the first argument is a function pointer
*/
extern void F_TRLAN(void(*op)(int*, int*, double*, int*, double*, int*),
	       int*, int*, int*, double*, double*, int*, double*, int*);

#ifdef __PGI
/* Portland Group Compilers need something special */
void pghpf_init(int *);
int __argc_save;
char **__argv_save;
static int zz = 0;
#endif

EXTERN_C DLLEXPORT int trlan(WolframLibraryData libData,
			mint Argc, MArgument *Args, MArgument Res) {

  mint I0;
  mint I1;
  int i, lwrk, nrow, mev, ipar[32];
  double tt, to, tp, tr;
  double eval[MEV], evec[NROW*MEV], wrk[MAXLAN*(MAXLAN+10)];


#ifdef __PGI
    __argc_save = argc;
    __argv_save = argv;
    pghpf_init(&zz);
#endif

    nrow = NROW; mev = MEV;
    lwrk = MAXLAN*(MAXLAN+10);
    ipar[0] = 0;
    ipar[1] = LOHI;
    ipar[2] = NED;
    ipar[3] = 0;     /* nec = 0 */
    ipar[4] = MAXLAN;
    ipar[5] = 5;     /* restarting scheme 1 */
    ipar[6] = 2000;  /* maximum number of MATVECs */
    ipar[7] = 0;     /* MPI_COMM */
    ipar[8] = -15;   /* verboseness */
    ipar[9] = 99;    /* Fortran IO unit number used to write log messages */
    ipar[10] = 1;    /* iguess -- 1 = use supply the starting vector */
    ipar[11] = -5;   /* checkpointing flag -- write about -5 times */
    ipar[12] = 98;   /* Fortran IO unit number to write checkpoint files */
    ipar[13] = 3*NROW;/* floating-point operations per MATVEC per PE */
    for (i=0; i<MEV; ++i) eval[i]=0.0;
    for (i=0; i<NROW; ++i) evec[i] = 1.0;
    wrk[0] = 1.4901E-8; /* relative tolerance on residual norms */

    /* call Fortran function TRLAN77 */
  printf("hi2\n");
    F_TRLAN(diag_op, ipar, &nrow, &mev, eval, evec, &nrow, wrk, &lwrk);
  printf("hi3\n");

    /* convert the time from clock ticks to seconds */
    tt = ((double) ipar[28]) * 1e-3; /* total time */
    tp = ((double) ipar[29]) * 1e-3; /* time in diag_op */
    to = ((double) ipar[30]) * 1e-3; /* time in orthogonalization */
    tr = ((double) ipar[31]) * 1e-3; /* time in restarting */
    /* printing the result -- trlan77 calls trl_print_info */
    if (my_pe == 0) {
	printf("MAXLAN:%10d, Restart:%10d,   NED:  -%7d,     NEC:%10d\n",
	       ipar[4], ipar[5], ipar[2], ipar[3]);
	printf("MATVEC:%10d,  Reorth:%10d, Nloop:%10d, Nlocked:%10d\n",
	       ipar[24], ipar[26], ipar[25], ipar[23]);
	printf("Ttotal:%10.6f,    T_op:%10.6f, Torth:%10.6f,  Tstart:%10.6f\n",
	       tt, tp, to, tr);
	for (i=0; i<ipar[3]; i++)
	    printf("E(%i) = %25.17g %16.4e\n", i+1, eval[i], wrk[i]);
    }

  I0 = MArgument_getInteger(Args[0]);
  I1 = I0 + 1;
  MArgument_setInteger(Res, I1);
  return LIBRARY_NO_ERROR;
}  /* end of trlan */

/*  The operator (matrix-vector multiplication routine)
    the matrix is a diagonal matrix with values (1, 4, 9, 16, 25, ....)
*/
void diag_op(int *nrow, int *ncol, double *xin, int *ldx, double *yout,
	     int *ldy) {
    int i, j, doff;
    double  *xtmp, *ytmp;

    for (j=0; j<*ncol; ++j) {
	xtmp = xin + *ldx * j;
	ytmp = yout + *ldy * j;
	for (i=1; i<=*nrow; ++i, ++xtmp, ++ytmp)
	    *ytmp = *xtmp * (i*i);
    }
} /* end of function diag_op */
