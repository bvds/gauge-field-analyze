(*
                     MINRES-QLP linear system solver

Import MATLAB code from https://web.stanford.edu/group/SOL/software/minresqlp/
The intention here is for the Mathematica code to match, as close as possible,
the original MATLAB code.

function [x, flag, iter, Miter, QLPiter, relres, relAres,
        Anorm, Acond, xnorm, Axnorm, resvec, Aresvec] =
    minresqlp(A, b, M, shift, rtol, maxit, maxxnorm, Acondlim,
        TranCond, show)

MINRESQLP: min-length solution to symmetric (possibly singular) Ax=b or
min||Ax-b||.

MINRESQLP(A,B) solves the system of linear equations A*X=B
or the least-squares problem min norm(B-A*X) if A is singular.
The N-by-N matrix A must be symmetric or Hermitian, but need not be
positive definite or nonsingular.  The rhs vector B must have length N.

MINRESQLP(AFUN,B) accepts a function handle AFUN instead of
the matrix A.  Y = AFUN(X) returns the matrix-vector product Y=A*X.
In all of the following syntaxes, A can be replaced by AFUN.

MINRESQLP(A,B,M) uses a matrix M as preconditioner.
M must be positive definite and symmetric or Hermitian.
It may be a function handle MFUN such that Y=MFUN(X) returns a
solution of M.Y=X.
If M is not supplied or M=None, a preconditioner is not applied.

MINRESQLP(A,B,M,SHIFT) solves (A - SHIFT*I)X = B, or the corresponding
least-squares problem if (A - SHIFT*I) is singular, where SHIFT is a
real scalar.
Default SHIFT = 0.

rToleraance specifies a stopping tolerance.
Default rTolerance = machine epsilon.

maxIterations specifies the maximum number of iterations.
Default maxIterations = 4*N.

maxXNorm, aConditionLimit, tranCondition are three parameters associated
with singular or ill-conditioned systems (A - SHIFT*I)*X = B.

maxXNorm is an upper bound on norm(X). Alternatively, maxXNorm may be
a function that acts on x, returning True if the bound is met.
Default maxXNorm = 10^7.

aConditionLimit is an upper bound on ACOND, an estimate of condition(A).
Default aConditionLimit = 10^15.

tranCondition is a real scalar >= 1.
If tranCondition>1,        a switch is made from MINRES iterations to
                      MINRES-QLP iterationsd when ACOND >= TRANCOND.
If tranCondition=1,        all iterations will be MINRES-QLP iterations.
If tranCondition=aConditionLimit, all iterations will be conventional
                   MINRES iterations (which are slightly cheaper).
Default trandCondition = 10^7.

printDetails specifies the printing level.
If printDetails>1,  an iteration log will be output.
If printDetails>0 or True,  a summary will be output.
If printDetails=0 or False, the log is suppressed.
Default printDetails=0.

If returnVectors is True, then the RESVEC and ARESVEC are returned.

Returns an array {X, FLAG, ITER, MITER, QLPITER, RELRES, RELARES,
ANORM, ACOND, XNORM, AXNORM, RESVEC, ARESVEC} where:

X is the solution

Convergence FLAG:
-1 (beta2=0)  B and X are eigenvectors of (A - SHIFT*I).
 0 (beta1=0)  B = 0.  The exact solution is X = 0.
 1 X solves the compatible (possibly singular) system (A - SHIFT*I)X = B
   to the desired tolerance:
      RELRES = RNORM / (ANORM*XNORM + NORM(B)) <= RTOL,
   where
           R = B - (A - SHIFT*I)X and RNORM = norm(R).
 2 X solves the incompatible (singular) system (A - SHIFT*I)X = B
   to the desired tolerance:
      RELARES = ARNORM / (ANORM * RNORM) <= RTOL,
   where
           AR = (A - SHIFT*I)R and ARNORM = NORM(AR).
 3 Same as 1 with RTOL = EPS.
 4 Same as 2 with RTOL = EPS.
 5 X converged to an eigenvector of (A - SHIFT*I).
 6 XNORM exceeded MAXXNORM.
 7 ACOND exceeded ACONDLIM.
 8 MAXIT iterations were performed before one of the previous
   conditions was satisfied.
 9 The system appears to be exactly singular.  XNORM does not
   yet exceed MAXXNORM, but would if further iterations were
   performed.

The number of iterations performed, with ITER = MITER + QLPITER.
 MITER   is the number of conventional MINRES iterations.
 QLPITER is the number of MINRES-QLP iterations.

The relative residuals for (A - SHIFT*I)X = B and the
associated least-squares problem.  RELRES and RELARES are
defined above in the description of FLAG.

Further norms:
    ANORM,  an estimate of the 2-norm of A-SHIFT*I.
    ACOND,  an estimate of COND(A-SHIFT*I,2).
    XNORM,  a recurred estimate of NORM(X).
    AXNORM, a recurred estimate of NORM((A-SHIFT*I)X)

Vectors:
    RESVEC,  a vector of estimates of NORM(R) at each iteration,
             including NORM(B) as the first entry.
    ARESVEC, a vector of estimates of NORM((A-SHIFT*I)R) at each
             iteration, including NORM((A-SHIFT*I)B) as the first entry.
    RESVEC and ARESVEC have length ITER+1.

See also BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, PCG, QMR, SYMMLQ,
TFQMR, CHOLINC, FUNCTION_HANDLE.
Also MINRES, SYMMLQ, LSQR, CGLS downloadable from
http://www.stanford.edu/group/SOL/software.html

REFERENCES:
Sou-Cheng T. Choi and Michael A. Saunders,
ALGORITHM: MINRES-QLP for Singular Symmetric and Hermitian Linear
Equations and Least-Squares Problems, to appear in ACM Transactions on
Mathematical Software.

Sou-Cheng T. Choi, Christopher C. Paige, and Michael A. Saunders,
MINRES-QLP: A Krylov Subspace Method for Indefinite or Singular Symmetric
Systems, SIAM Journal on Scientific Computing, Vol. 33, No. 4, August
2011, pp. 1810--1836.

Sou-Cheng T. Choi's PhD Dissertation, Stanford University, 2006:
http://www.stanford.edu/group/SOL/dissertations.html

 CURRENT / FUTURE RELEASES of minresqlp:
Version 2: 
   http://code.google.com/p/minres-qlp/
   http://www.mathworks.com/matlabcentral/fileexchange
Version 1: 
   http://code.google.com/p/minres-qlp/
   http://www.stanford.edu/group/SOL/download.html
Other implementations in Fortran 90/95, Python:
   http://code.google.com/p/minres-qlp/


MODIFICATION HISTORY:
28 Jun 2013: Second version for MATLAB Central.
27 Jun 2013: (1) Fixed iteration log:
             (a) The heading came out ok every 20 lines initially, but
             stopped after itn 80.
             (b) Subsequent output was for itns 109, 119, ... rather
             than 110, 120,...
             (2) Introduced local variable "likeLS", which is true if
             Ax = b looks more like a least-squares system.
             (3) Fixed an error in minresxxxM(). Moved two lines that
             involve R into the else statement. Added code for handling
             preconditioner as a function handle.
             (4) Added debug statements.
28 Jul 2012: (1) Fixed a bug in pnorm_0 used in Anorm recurrence
             relation so that {Anorm(k)} are monotonic increasing
             underestimates of ||A||_2.
             (2) Fixed SymOrtho to ensure the 2x2 Hermitian reflectors
             are orthonormal.
20 Jul 2012: (1) Changed default RTOL to 1e-15 and default MAXIT to 4N.
             (2) Initalized a few local variables to zeros.
12 May 2006: Created MINRESQLPs.m from research file minresqlp35.m.
19 Apr 2010: Help formatted like Matlab routines. Added output Aresvec.
25 Apr 2010: Eliminated unreferenced variables. Shortened
             resvec and Aresvec from maxit to (iter+1).
             The final rnorm, Arnorm, xnorm are computed directly.
26 Apr 2010: Special tests for beta2=0 and/or alfa1=0 when no
             preconditioning.
02 May 2010: First version for
             http://www.stanford.edu/group/SOL/software.html.

AUTHORS: Sou-Cheng (Terrya) Choi, CI, University of Chicago
       Michael Saunders, SOL, Stanford University

COPYRIGHT NOTICE:

This is Copyrighted Material. The software is COPYRIGHTED by the
original authors.

COPYLEFT NOTICE:

Permission is granted to make and distribute verbatim copies of this
file, provided that the entire file is copied **together** as a
**unit**.

The purpose of this permission notice is to allow you to make copies
of the software and distribute them to others, for free or for a fee,
subject to the constraint that you maintain the entire file here
as a unit.  This enables people you give the software to be aware
of its origins, to ask questions of us by e-mail, to request
improvements, obtain later releases, and so forth.

If you seek permission to copy and distribute translations of this
software into another language, please e-mail a specific request to
saunders@stanford.edu and scchoi@stanford.edu.

If you seek permission to excerpt a **part** of the software library,
for example to appear in a scientific publication, please e-mail a
specific request to saunders@stanford.edu and scchoi@stanford.edu.

COMMENTS?

Email sctchoi@uchicago.edu and saunders@stanford.edu
*)

minresqlp::usage = "MINRES-QLP solver for symmetric indefinte linear systems.  The matrix can be expressed as a pure function that acts on a vector.";
minresqlp::indefinite = "`1` appears to be indefinite.";
Options[minresqlp] = {rTolerance -> $MachineEpsilon,
   maxIterations -> Automatic, maxXNorm -> 10.0^7, aConditionLimit -> 10.0^15,
   tranCondition -> 10.0^7, printDetails -> False, returnVectors->True};
minresqlp[A_?MatrixQ, rest__] := minresqlp[(A.#)&, rest];
minresqlp[A_, b_, M_?MatrixQ, rest___] :=
  minresqlp[A, b, LinearSolve[M], rest];
minresqlp[A_Function, b_, M:(None|_Function|_LinearSolveFunction):None,
	  shift:_?NumberQ:0, OptionsPattern[]] :=
Module[{
    (* MATLAB constants *)
    eps = $MachineEpsilon,
    realmax = $MaxMachineNumber,
    realmin = $MinMachineNumber,
    (* Mimic MATLAB function *)
    zeros = Function[Array[0.0&, {##}]],

    show = Replace[OptionValue[printDetails], {False -> 0, True -> 1}],
    rtol = OptionValue[rTolerance], maxxnorm = OptionValue[maxXNorm],
    Acondlim = OptionValue[aConditionLimit],
    TranCond = OptionValue[tranCondition],
    maxit = OptionValue[maxIterations],

    debug = False,
    n = Length[b],
    resvec = None, Aresvec = None, r2, r3, beta1},
If[maxit === Automatic, maxit = 4 n];
If[OptionValue[returnVectors],      
   resvec = zeros[maxit + 1];
   Aresvec = zeros[maxit + 1]];

(* Set up {beta1,p,v} for the first Lanczos vector v1. *)
r2 = N[b];        (* r2=b *)
r3 = r2;          (* r3=b *)
beta1 = Norm[r2]; (* beta1=norm(b) *)
If[M =!= None,
   r3 = M[r2]; (* M*r3=b *)
   beta1 = Conjugate[r3].r2; (* beta1=b'*inv(M)*b *)
   If[beta1 < 0,
      Message[minresqlp::indefinite, "M"];
      beta1 = $Failed,
      beta1 = Sqrt[beta1]]];

(* Initialize other quantities. *)
Block[{
flag0 = -2, flag,
iter = 0, QLPiter = 0,
lines = 1, headlines = 20,
beta = 0, tau = 0, taul = 0, phi = beta1,
betan = beta1, gmin = 0, cs = -1, sn = 0,
cr1 = -1, sr1 = 0, cr2 = -1, sr2 = 0,
dltan = 0, eplnn = 0, gama = 0, gamal = 0,
gamal2, eta = 0, etal = 0, etal2 = 0,
vepln = 0, veplnl = 0, veplnl2 = 0, ul3 = 0,
ul2 = 0, ul = 0, u, rnorm,
xnorm = 0, xl2norm = 0, Axnorm = 0,
Anorm = 0, Acond = 1,
gamaQLP = 0, gamalQLP = 0, veplnQLP = 0, gminl = 0,
uQLP = 0, ulQLP = 0,
relres,
relresl = 0,
relAresl = 0,
v, x, xl2, w, wl, wl2, r1,
betal, alfa, pnorm, dlta, gbar, xnorml, taul2, epln,
likeLS, Anorml, Acondl, rnorml, Arnorml, dltaQLP, gamaTmp, gamalTmp},
flag = flag0; rnorm = betan;
relres = rnorm/(beta1 + 10.0^-50); (* Safeguard for beta1=0 *)
x = zeros[n];
xl2 = x;
w = x;
wl = x;
r1 = x;

If[resvec =!= None, resvec[[1]] = Re[beta1]];

(* print header if show *)
Block[{first = "Enter minresqlp.  ",
last = "Exit minresqlp.  ",
msg = {
        " beta2=0.  b and x are eigenvectors                   ", (* -1 *)
        " beta1=0.  The exact solution is  x = 0               ", (* 0 *)
        " A solution to Ax = b found, given rtol                 ", (* 1 *)
        " Min-length solution for singular LS problem, given rtol", (* 2 *)
        " A solution to Ax = b found, given eps                  ", (* 3 *)
	" Min-length solution for singular LS problem, given eps ", (* 4 *)
        " x has converged to an eigenvector                      ", (* 5 *)
	" xnorm has exceeded maxxnorm                            ", (* 6 *)
        " Acond has exceeded Acondlim                            ", (* 7 *)
	" The iteration limit was reached                        ", (* 8 *)
        " Least-squares problem but no converged solution yet    "}, (* 9 *)
head = {"iter", "rnorm", "Arnorm", "Compatible", "LS",
	"Anorm","Acond","xnorm"}},
If[show > 1,
   Print[first];
   Print["Min-length solution of symmetric (A-sI)x = b or min ||(A-sI)x - b||"];
   Print["n=", n, "   ||b||=", beta1, "   shift=", shift,
	 "   rtol=",rtol];
   Print["maxit=", maxit, "   maxxnorm=", If[NumberQ[maxxnorm],maxxnorm,
					   "function"],
	 "   Acondlim=", Acondlim, "   TranCond=", TranCond];
   Print[" ", head]];

If[beta1 == 0, flag = 0]; (* b=0 => x=0. We will skip the main loop. *)

(* Main iteration *)
While[flag == flag0 && iter < maxit,

    (* Lanczos *)
    iter = iter + 1;
    betal = beta; beta = Re[betan];
    v = r3*(1/beta); r3 = A[v];
    If[shift != 0, r3 = r3 - shift*v];
    If[iter > 1, r3 = r3 - (beta/betal)*r1];

    alfa = Re[Conjugate[r3].v]; (* Allow for Hermitian A. Must get real alfa here. *)
    r3 = r3 - (alfa/beta)*r2; r1 = r2; r2 = r3;

    If[M === None,
       betan = Norm[r3];
       If[iter == 1, (* Something special can happen *)
          If[betan == 0, (* beta2=0 *)
             If[alfa == 0, (* alfa1=0 *)
		flag = 0; (* Ab = 0 and x = 0  ("A" = (A - shift*I)) *)
		Break[],
		flag = -1; (* Ab = alfa1 b, x = b/alfa1, an eigenvector *)
		x = N[b]/alfa;
		Break[]]]],
       r3 = M[r2]; betan = Conjugate[r2].r3;
       If[betan > 0,
	  betan = Sqrt[betan],
          Message[minresqlp::indefinite, "M"];
          betan = $Failed]];
    If[iter <= 2,
       pnorm = Norm[{alfa, betan}],
       pnorm = Norm[{betal, alfa, betan}]];

    If[debug,
       Print["Lanczos iteration ", iter, ":"];
       Print["  v_", iter, "     = ", Take[v, Min[n, 5]]];
       Print["  r1_", iter, "    = ", Take[r1, Min[n, 5]]];
       Print["  r2_", iter, "    = ", Take[r2, Min[n, 5]]];
       Print["  r3_", iter, "    = ", Take[r3, Min[n, 5]]];
       Print["  alpha_", iter, " = ", alfa, ", beta_", iter, " = ",
             beta, ", beta_", iter + 1, " = ", betan, " pnorm_", iter,
             " = ", pnorm]];

    (* Apply previous left reflection Q_{k-1} *)
    Block[{dbar = dltan},
    dlta = cs*dbar + sn*alfa; epln = eplnn;
    gbar = sn*dbar - cs*alfa; eplnn = sn*betan;
    dltan = -cs*betan; dltaQLP = dlta;

    If[debug,
       Print["Apply previous left reflection Q_{", iter-1, ",", iter, "}:"];
       Print["  c_", iter-1, "     = ", cs,", s_", iter-1,"    = ",sn];
       Print["  dlta_", iter, " = ", dlta, ", gbar_", iter, " = ", gbar];
       Print["  epln_", iter+1, " = ", eplnn, ", dbar_", iter+1," = ",dltan]]];

    (* Compute the current left reflection Q_k *)
    gamal2 = gamal; gamal = gama;
    {cs, sn, gama} = SymOrtho[gbar, betan]; gamaTmp = gama;
    taul2 = taul; taul = tau; tau = cs*phi;
    Axnorm = Norm[{Axnorm, tau}]; phi = sn*phi;

    If[debug,
       Print["Compute the current left reflection Q_{",iter,",",iter+1,"}:"];
       Print["  c_",iter, "     = ",cs, ", s_",iter,"    = ",sn];
       Print["  tau_",iter,"   = ",tau,", phi_",iter,"  = ",phi];
       Print["  gama_", iter, " = ", gama]];

    (* Apply the previous right reflection P{k-2,k} *)
    Block[{dltaTmp},
    If[iter > 2,
       veplnl2 = veplnl; etal2 = etal; etal = eta;
       dltaTmp = sr2*vepln - cr2*dlta;
       veplnl = cr2*vepln + sr2*dlta;
       dlta = dltaTmp; eta = sr2*gama; gama = -cr2*gama;

       If[debug,
          Print["Apply the previous right reflection P_{",iter-2,",",iter,"}:"];
          Print["  cr2_",iter, "   = ",cr2,", sr2_",iter, "    = ",sr2];
          Print["  gama_",iter-2, " = ", gamal2, ", gama_",iter-1,
          "  = ",gamal,", gama_",iter," = ",gama];
          Print["  dlta_",iter, " = ", dlta, ", vepln_",iter-1,
		" = ",veplnl, ", eta_",iter,"   = ", eta]]]];

    (* Compute the current right reflection P{k-1,k}, P_12, P_23, ... *)
    If[iter > 1,
       {cr1, sr1, gamal} = SymOrtho[gamal, dlta];
       vepln = sr1*gama;
       gama = -cr1*gama;

       If[debug,
          Print["Compute the second current right reflections P_{",iter-1,",",iter,"}:"];
          Print["  cr1_",iter, "   = ", cr1, ", sr1_",iter,"   = " sr1];
          Print["  gama_",iter-1, " = ", gamal, ", gama_",iter,
		" = ", gama, ", vepln_",iter, " = ", vepln]]];

    (* Update xnorm *)
    Block[{ul4 = ul3,xnormTmp},
    xnorml = xnorm; ul3 = ul2;
    If[iter > 2,
       ul2 = (taul2 - etal2*ul4 - veplnl2*ul3)/gamal2];
    If[iter > 1,
       ul = (taul - etal*ul3 - veplnl*ul2)/gamal];
    xnormTmp = Norm[{xl2norm, ul2, ul}];
    likeLS = (relresl >= relAresl);
    If[Abs[gama] > realmin && xnormTmp < maxxnorm,
       u = (tau - eta*ul2 - vepln*ul)/gama;
       If[Norm[{xnormTmp, u}] > maxxnorm && likeLS,
	  u = 0; flag = 6],
       u = 0; flag = 9];
    xl2norm = Norm[{xl2norm, ul2}];
    xnorm = Norm[{xl2norm, ul, u}]];

    (* Update w. Update x except if it will become too big *)
    If[Acond < TranCond && flag == flag0 && QLPiter == 0,
       (* MINRES updates *)
       wl2 = wl; wl = w;
       w = (v - epln*wl2 - dltaQLP*wl)*(1/gamaTmp);
       If[xnorm < maxxnorm,
	  x = x + tau*w,
	  flag = 6],

       (* MINRES-QLP updates *)
       QLPiter = QLPiter + 1;
       If[QLPiter == 1,
	  xl2 = zeros[n];
          If[iter > 1, (* construct w_{k-2}, w_{k-1} *)
             If[iter > 2,
		wl = gamalQLP*wl + veplnQLP*w]; (* w_{k-2} *)
             w = gamaQLP*w; xl2 = x - wl*ulQLP - w*uQLP]];
       Which[iter == 1,
	     wl2 = wl; wl = v*sr1; w = -v*cr1,
	     iter == 2,
             wl2 = wl;
             wl = w*cr1 + v*sr1;
             w = w*sr1 - v*cr1,
             True,
	     wl2 = wl; wl = w; w = wl2*sr2 - v*cr2;
             wl2 = wl2*cr2 + v*sr2; v = wl*cr1 + w*sr1;
             w = wl*sr1 - w*cr1; wl = v];
       xl2 = xl2 + wl2*ul2;
       x = xl2 + wl*ul + w*u];

    If[debug,
       Print["Update w:"];
       Print["  w_",iter, "     = ", Take[wl, Min[n, 5]]];
       Print["  w_",iter, "     = ", Take[w, Min[n, 5]]];
       Print["Update u, x and xnorm:"];
       Print["  u_",iter - 2, "     = ", ul2, ", u_",iter - 1,
             "     = ", ul, ", u_",iter, "     = ", u];
       Print["  x_",iter, "     = ", Take[x, Min[n, 5]]];
       Print["  ||x_",iter, "|| = ", xnorm]];

    (* Compute the next right reflection P{k-1,k+1} *)
    gamalTmp = gamal;
    {cr2, sr2, gamal} = SymOrtho[gamal, eplnn];

    (* Store quantities for transfering from MINRES to MINRES-QLP *)
    gamalQLP = gamalTmp; veplnQLP = vepln; gamaQLP = gama;
    ulQLP = ul; uQLP = u;

    (* Estimate various norms *)
    Block[{gminl2, rootl, absGama},
    absGama = Abs[gama]; Anorml = Anorm;
    Anorm = Max[{Anorm, gamal, absGama, pnorm}];
    Which[iter == 1,
	  gmin = gama; gminl = gmin,
	  iter > 1,
          gminl2 = gminl; gminl = gmin; gmin = Min[{gminl2, gamal, absGama}]];
    Acondl = Acond; Acond = Anorm/gmin;
    rnorml = rnorm; relresl = relres;
    If[flag != 9,
       rnorm = phi];
    relres = rnorm/(Anorm*xnorm + beta1);
    rootl = Norm[{gbar, dltan}];
    Arnorml = rnorml*rootl;
    relAresl = rootl/Anorm];

    (* See if any of the stopping criteria are satisfied. *)
    Block[{epsx = Anorm*xnorm*eps, t1, t2},
	  If[flag == flag0 || flag == 9,
	     t1 = 1 + relres;
	     t2 = 1 + relAresl;
             If[iter >= maxit, flag = 8]; (* Too many itns *)
             If[Acond >= Acondlim, flag = 7]; (* Huge Acond *)
             If[xnorm >= maxxnorm, flag = 6]; (* xnorm exceeded its limit *)
             If[epsx >= beta1, flag = 5]; (* x is an eigenvector *)
             If[t2 <= 1, flag = 4]; (* Accurate LS solution *)
             If[t1 <= 1, flag = 3]; (* Accurate Ax=b solution *)
             If[relAresl <= rtol, flag = 2]; (* Good enough LS solution *)
             If[relres <= rtol, flag = 1]; (* Good enough Ax=b solution *)

             If[debug,
		Print["Update other norms:"];
		Print["  gmin_",iter, "   = ", gmin];
		Print["  pnorm_",iter, "  = ", pnorm];
		Print["  rnorm_",iter, "  = ", rnorm];
		Print["  Arnorm_",iter - 1, " = ", Arnorml];
		Print["  Acond_",iter, "  = ", Acond]]]];

(*
The "disable" option allowed iterations to continue until xnorm
became large and x was effectively a nullvector.
We know that r will become a nullvector much sooner,
so we now disable the disable option.  :-)

If[disable && (iter<maxit),
    flag=0;
    If[Axnorm<rtol*Anorm*xnorm,
        flag=10]];
*)

    If[flag == 2 || flag == 4 || (flag == 6 && likeLS) || flag == 7, (* Possibly singular *)
       iter = iter - 1;
       Acond = Acondl; rnorm = rnorml; relres = relresl,
       If[resvec =!= None,
	  resvec[[iter + 1]] = rnorm;
          Aresvec[[iter]] = Arnorml];
       If[show > 1 && Mod[iter-1,lines] == 0,
          Which[iter == 101,
		lines = 10; headlines = 20*lines,
		iter == 1001,
		lines = 100; headlines = 20*lines];
          Print[If[QLPiter == 1, "P", " "],
		{iter-1, rnorml, Arnorml,
		 relresl, relAresl, Anorml, Acondl, xnorml}];
          If[iter > 1 && Mod[iter, headlines] == 1,
	     Print[head]]]]];

(* We have exited the main loop. *)
Block[{Arnorm, relAres,
       start = If[QLPiter == 1, "P", " "],
       Miter = iter - QLPiter},

(* Compute final quantities directly. *)
r1 = b - A[x] + shift*x; (* r1 is workspace for residual vector *)
rnorm = Norm[r1];
Arnorm = Norm[A[r1] - shift*r1];
xnorm = Norm[x];
relres = rnorm/(Anorm*xnorm + beta1);
relAres = 0;
If[rnorm > realmin, relAres = Arnorm/(Anorm*rnorm)];
If[Aresvec =!= None,
   Aresvec[[iter + 1]] = Arnorm;
   Aresvec = Take[Aresvec,iter + 1]];
If[resvec =!= None,
   resvec = Take[resvec,iter + 1]];
If[show > 1,
   If[rnorm > realmin,
      Print[start, {iter, rnorm, Arnorm, relres, relAres, Anorm,
		    Acond, xnorm}],
      Print[start, {iter, rnorm, Arnorm, relres, 0, Anorm, Acond,
		    xnorm}]]];
If[show > 0,
   Print[last, " flag=", flag, "   ", msg[[flag + 2]]];
   Print[last, " iter=", iter, "   (MINRES ", Miter, ", MINRES-QLP ", QLPiter, ")"];
   Print[last, " rnorm=", rnorm, "   rnorm direct=", Norm[r1]];
   Print[last, "                   Arnorm direct=", Arnorm];
   Print[last, " xnorm=", xnorm, "   xnorm direct=", Norm[x]];
   Print[last, " Anorm=", Anorm, "   Acond=", Acond]];

(* Return values *)
{x, flag, iter, Miter, QLPiter, relres, relAres, Anorm, Acond, xnorm, Axnorm,
 resvec, Aresvec}]]]];


(* 
  SymOrtho: Stable Symmetric Householder reflection

INPUTS:
 a      first element of a two-vector  (a, b)     
 b      second element of a two-vector (a, b)

OUTPUTS:
 c      cosine(theta), where theta is the implicit angle of rotation
        (counter-clockwise) in a plane-rotation
 s      sine(theta)
 r      two-norm of (a, b)
 
DESCRIPTION:
  Stable symmetric Householder reflection that gives c and s such that
     ( c  s )(a) = (d),
     ( s -c )(b) = (0) 
  where d = two-norm of vector (a, b),
     c = a / sqrt(a^2 + b^2) = a / d,
     s = b / sqrt(a^2 + b^2) = b / d.
  The implementation guards against overlow in computing sqrt(a^2 + b^2).

EXAMPLE:
   description

SEE ALSO:
  TESTSYMGIVENS.m,
  PLANEROT (MATLAB's function) --- 4 divisions while 2 would be enough,
  though not too time-consuming on modern machines
 
REFERENCES:
 Algorithm 4.9, stable *unsymmetric* Givens rotations in
  Golub and van Loan's book Matrix Computations, 3rd edition.

MODIFICATION HISTORY:
 10/06/2004: Replace d = norm([a,b]) by
                     d = a/c if |b| < |a| or b/s otherwise.
 10/07/2004: First two cases (either b or a == 0) rewritten to make sure
             (1) d >= 0
             (2) if [a,b] = 0, then c = 1 and s = 0 but not c = s = 0.
 09/27/2011: Change filename from SYMGIVENS2 to SYMORTHO.
 01/16/2012: Change file from SYMORTHO to SYMREFL.
 

KNOWN BUGS:
  MM/DD/2004: description

AUTHORS: Sou-Cheng (Terrya) Choi, CI, University of Chicago
       Michael Saunders, SOL, Stanford University

CREATION DATE: 09/28/2004
*)

SymOrtho[a_, b_] := Block[{
   
absa = Abs[a],
absb = Abs[b],
signa = Sign[a],
signb = Sign[b],
t,c,s,r},

If[ Im[a]==0 && Im[b]==0, 
    (* Both a and b are real numbers *)

    (* Special cases: a or b is 0 *)
    Which[b == 0,
	  If[a == 0,
	     c = 1.0,
	     c = signa]; (* NOTE: Sign(0) = 0 in Mathematica *)
	  s = 0;
	  r = absa,

	  a == 0,
	  c = 0;
	  s = signb;
	  r = absb,

	  (* Both a and b are non-zero *)
	  absb > absa,
	  t = a/b;
	  s = signb / Sqrt[1 + t^2];
	  c = s*t;
	  r = b/s, (* computationally better than d = a / c since |c| <= |s| *)

	  True,
	  t = b/a;
	  c = signa / Sqrt[1 + t^2];
	  s = c*t;
	  r = a/c], (* computationally better than d = b / s since |s| <= |c| *)

    (* a and/or b are complex numbers *)
    (* Special cases: a or b is 0 *)
    Which[ b == 0,
	   c = 1;
	   s = 0;
	   r = a,

	   a == 0,
	   c = 0;
	   s = 1;
	   r = b,

	   (* Both a and b are non-zero *)
	   absb > absa,
	   t = absa/absb;
	   c = 1/Sqrt[1+t^2]; (* temporary *)
	   s = c*Conjugate[signb/signa];
	   c = c*t;
	   r = b/Conjugate[s],

	   True,
	   t = absb/absa;
	   c = 1/Sqrt[1+t^2];
	   s = c*t*Conjugate[signb/signa];
	   r = a/c]];

{c, s, r}];
