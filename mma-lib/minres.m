(*
                     MINRES linear system solver
 
Import MATLAB code from https://web.stanford.edu/group/SOL/software/minres/
The intention here is for the Mathematica code to match, as close as possible,
the original MATLAB code.
 
function [ x, istop, itn, rnorm, Arnorm, Anorm, Acond, ynorm resvec] = ...
           minres( A, b, M, shift, printDetails, check, maxIterations,
                   stepMonitor, rTolerance, localSize )

minres solves the n x n system of linear equations Ax = b
or the n x n least squares problem           min ||Ax - b||_2^2,
where A is a symmetric matrix (possibly indefinite or singular)
and b is a given vector.  The dimension n is defined by length(b).

INPUT:

"A" may be a dense or sparse matrix (preferably sparse!)
or a function handle such that y = A(x) returns the product
y = A*x for any given n-vector x.

If "M" = None or not supplied, preconditioning is not used.  Otherwise,
"M" defines a positive-definite preconditioner M = C*C^T.
"M" may be a dense or sparse matrix (preferably sparse!)
or a function handle such that y = M(x) solves the system
My = x given any n-vector x.

If shift != 0, minres really solves (A - shift*I)x = b
(or the corresponding least-squares problem if shift is an
eigenvalue of A).

When M = C*C^T exists, minres implicitly solves the system

           P(A - shift*I)P^T xbar = Pb,
   i.e.               Abar xbar = bbar,
   where                      P = inv(C),
                           Abar = P(A - shift*I)P^T,
                           bbar = Pb,

and returns the solution      x = P^T xbar.
The associated residual is rbar = bbar - Abar xbar
                                = P(b - (A - shift*I)x)
                                = Pr.

OUTPUT:

x      is the final estimate of the required solution
       after k iterations, where k is return in itn.
istop  is a value from [-1:9] to indicate the reason for termination.
       The reason is summarized in msg[istop+2] below.
itn    gives the final value of k (the iteration number).
rnorm  estimates norm(r_k)  or norm(rbar_k) if M exists.
Arnorm estimates norm(Ar_{k-1}) or norm(Abar rbar_{k-1}) if M exists.
       NOTE THAT Arnorm LAGS AN ITERATION BEHIND rnorm.

resvec returns a vector of the residual at each iteration

Code author: Michael Saunders, SOL and ICME, Stanford University
Contributors:Chris Paige, School of Computer Science, McGill University
             Sou-Cheng Choi, ICME, Stanford University
 

Known bugs:
  1. As Jeff Kline pointed out, Arnorm = ||A r_{k-1}|| lags behind
     rnorm = ||r_k||.  On singular systems, this means that a good
     least-squares solution exists before Arnorm is small enough
     to recognize it.  The solution x_{k-1} gets updated to x_k
     (possibly a very large solution) before Arnorm shuts things
     down the next iteration.  It would be better to keep x_{k-1}.
*)

minres::usage = 
  "MINRES solver for symmetric indefinte linear systems.  The matrix can be expressed as a pure function that acts on a vector.";
Options[minres] := {printDetails -> 0, check -> True,
    maxIterations -> Automatic, rTolerance -> $MachineEpsilon,
    localSize -> 0, stepMonitor -> None};
minres[A_?MatrixQ, rest__] := minres[(A.#) &, rest];
minres[A_, b_, M_?MatrixQ, rest___] := 
  minres[A, b, LinearSolve[M], rest];
minres[A_Function, b_, 
  M:(None|_Function|_LinearSolveFunction):None, 
  shift:_?NumberQ:0.0, OptionsPattern[]] := 

 (* Initialize *)
 Module[{(* MATLAB constants *) eps = $MachineEpsilon, 
   realmax = $MaxMachineNumber, (* Mimic MATLAB function *)
   zeros = Function[Array[0.0 &, {##}]],
   localOrtho, localVEnqueue, localVOrtho,
   last = "Exit minres.  ", 
   msg = {"beta2=0. If M=I, b and x are eigenvectors", 
     "beta1=0. The exact solution is x=0", 
     "A solution to Ax=b was found, given rtol", 
     "A least-squares solution was found, given rtol" , 
     "Reasonable accuracy achieved, given machine epsilon", 
     "x has converged to an eigenvector" , 
     "acond has exceeded 0.1/machine epsilon", 
     "The iteration limit was reached", 
     "A does not define a symmetric matrix", 
     "M does not define a symmetric matrix", 
	  "M does not define a pos-def preconditioner"},
    n = Length[b],	 
   itnlim = OptionValue[maxIterations], rtol = OptionValue[rTolerance], 
   show = Replace[OptionValue[printDetails], {False -> 0, True -> 1}],
   istop = 0, itn = 0, Anorm = 0, Acond = 0,
   rnorm = 0, ynorm = 0, done = False, x, resvec, y, r1, beta1, Arnorm,
   tinit = SessionTime[], ta = 0, tm = 0, tortho = 0,
   addTime = Function[{timer, expr},
     Block[{t1 = SessionTime[], result = expr}, 
	   timer += SessionTime[] - t1; result], {HoldAll, SequenceHold}]},
  If[itnlim === Automatic, itnlim = 5 n]; (* Not in MATLAB code *)
  If[show > 1,
    Print["minres.m SOL, Stanford University Version of 2015"];
    Print["Solution of symmetric Ax=b or (A-shift*I)x=b"];
    Print["n=", n, " shift=", shift];
    Print["itnlim=", itnlim, " rtol=", rtol]];
  x = zeros[n];
  resvec = zeros[itnlim];

  (* Initialization for local reorthogonalization *)
  {localOrtho, localVEnqueue, localVOrtho} = 
   makeOrtho[n, OptionValue[localSize]];

  (* Set up y and v for the first Lanczos vector v1.
     y=beta1 P^T v1, where P=C**(-1).
     v is really P^T v1. *)
  y = N[b];
  r1 = N[b]; (* initial guess x=0 initial residual *)
  If[M =!= None, addTime[tm, y = M[b]]];
  beta1 = Conjugate[b].y;

  (* Test for an indefinite preconditioner.
  If b=0 exactly, stop with x=0. *)

  If[beta1 < 0, istop = 9; show = 1; done = True];
  If[beta1 == 0, show = 1; done = True];
  If[beta1 > 0,
   beta1 = Sqrt[beta1]; (* Normalize y to get v1 later.*)

   (* See if M is symmetric. *)
   If[OptionValue[check] && M =!= None, 
    Block[{r2, s, t, z, epsa},
     r2   = M[y];
     s    = y.y;
     t    = r1.r2;
     z    = Abs[s - t];
     epsa = (s + eps)*eps^(1/3);
     If[z > epsa, istop = 8; show = 1; done = True]]];

   (* See if A is symmetric. *)
   If[OptionValue[check], 
    Block[{r2, w, s, t, z, epsa},
     w    = A[y];
     r2   = A[w];
     s    = w.w;
     t    = y.r2;
     z    = Abs[s - t];
     epsa = (s + eps)*eps^(1/3);
     If[ z > epsa, istop = 7; done = True; show = 1]]]];

  (* Initialize other quantities. *)
  Block[{oldb = 0, beta = beta1, dbar = 0, epsln = 0, qrnorm = beta1, 
    phibar = beta1, rhs1 = beta1, rhs2 = 0, tnorm2 = 0, gmax = 0, 
    gmin = realmax, cs = -1, sn = 0,
    w = zeros[n],
    w2 = zeros[n], 
    r2 = r1, alfa}, 
   If[show > 1,
      Print[{"Itn", "x(1)", "Compatible", "LS", "norm(A)", "cond(A)",
         "gbar/|A|" (* Check gbar *)}]];

   (* Main iteration loop. *)
   If[! done, (* k=itn=1 first time through *)
    While[itn < itnlim, (* max num of iter *)
     itn = itn + 1;

     (* Obtain quantities for the next Lanczos vector vk+1, k=1, 2, ...

        The general iteration is similar to the case k=1 with v0=0:
         p1=Operator*v1-beta1*v0,
         alpha1=v1^Tp1,q2=p1-alpha1*v1,
         beta2^2=q2^Tq2,
         v2=(1/beta2) q2.
     Again, y=betak P vk, where P=C**(-1).
     ....more description needed. *)
     s = 1/beta; (* Normalize previous vector (in y).*)
     v = s*y; (* v=vk if P=I *)

     (* if localOrtho turned on store old v for local reorthogonaliztion of new v *)
     If[localOrtho,
	addTime[tortho, localVEnqueue[v]]];
     addTime[ta, y = A[v] - shift*v]; (* shift is 0 otherwise solving A-shift*I *)
     If[itn >= 2,
	y = y - (beta/oldb)*r1]; (* normalization is the division r1 by oldb *)
     alfa = v.y; (* alphak *)
     y = (-alfa/beta)*r2 + y; (* normalization of r2/beta=v *)
     If[localOrtho,
	(* v will be normalized through y later- this is explicit
	 orthogonalizing versus the previous localSize lanczos vectors *)
	addTime[tortho, y = localVOrtho[y]]];
     r1 = r2; (* r1 is unnormalized vold *)
     r2 = y; (* r2 is unnormalized v *)
     If[M =!= None, addTime[tm, y = M[r2]]];
     oldb = beta; (* oldb=betak *)
     beta = r2.y; (* beta=betak+1^2 *)
     If[beta < 0, istop = 9; Break[]];
     beta = Sqrt[beta];
     tnorm2 = tnorm2 + alfa^2 + oldb^2 + beta^2;

     If[itn == 1,  (* Initialize a few things. *)
      If[beta/beta1 <= 10*eps, (* beta2=0 or~0. *)
	 istop = -1 (* Terminate later. *)]];

     (* Apply previous rotation Qk-1 to get
        [deltak epslnk+1]=[cs sn][dbark 0]
        [gbar k dbar k+1][sn-cs][alfak betak+1]. *)
     Block[{oldeps, delta, gbar, root, gamma, phi, epsx, epsr}, 
      oldeps = epsln;
      delta = cs*dbar + sn*alfa; (* delta1=0 deltak *)
      gbar  = sn*dbar - cs*alfa; (* gbar 1=alfa1 gbar k *)
      epsln = sn*beta; (* epsln2=0 epslnk+1 *)
      dbar  = -cs*beta; (* dbar 2=beta2 dbar k+1 *)
      root = Norm[{gbar, dbar}];
      Arnorm = phibar*root; (* ||Ar{k-1}|| *)

      (* Compute the next plane rotation Qk *)
      gamma = Norm[{gbar, beta}]; (* gammak *)
      gamma = Max[gamma, eps];
      cs = gbar/gamma; (* ck *)
      sn = beta/gamma; (* sk *)
      phi = cs*phibar; (* phik *)
      phibar = sn*phibar; (* phibark+1 *)

      (* Update x. *)
      Block[{denom, w1},
      denom = 1/gamma;
      w1 = w2;
      w2 = w;
      w = (v - oldeps*w1 - delta*w2)*denom];
      x = x + phi*w;

      If[OptionValue[stepMonitor] =!= None,
	 OptionValue[stepMonitor][itn, x, w, v]];

      (* Go round again. *)
      gmax = Max[gmax, gamma];
      gmin = Min[gmin, gamma];
      Block[{z = rhs1/gamma},
       rhs1 = rhs2 - delta*z;
       rhs2 = -epsln*z];

      (* Estimate various norms. *)
      Block[{epsa, diag},
       Anorm = Sqrt[tnorm2];
       ynorm = Norm[x];
       epsa = Anorm*eps;
       epsx = Anorm*ynorm*eps;
       epsr = Anorm*ynorm*rtol;
       diag = gbar;
       If[diag == 0, diag = epsa]];

      qrnorm = phibar;
      rnorm = qrnorm;
      resvec[[itn]] = rnorm;
      Block[{test1, test2},
	test1 = rnorm/(Anorm*ynorm); (* ||r|| /(||A|| ||x||) *)
	test2 = root/Anorm; (* ||Ar{k-1}|| /(||A|| ||r_{k-1}||) *)

       (* Estimate cond(A).
	  In this version we look at the diagonals of R in the
	  factorization of the lower Hessenberg matrix, Q*H=R,
          where H is the tridiagonal matrix from Lanczos with one
	  extra row, beta(k+1) e_k^T. *)
       Acond = gmax/gmin;

       (* See if any of the stopping criteria are satisfied. 
          In rare cases, istop is already -1 from above (Abar=const*I). *)
       If[istop == 0, 
        Block[{t1, t2},
	 t1 = 1 + test1; (* These tests work if rtol<eps *)
	 t2 = 1 + test2;
	 If[t2 <= 1, istop = 2];
         If[t1 <= 1, istop = 1];
         If[itn >= itnlim, istop = 6];
         If[Acond >= 0.1/eps, istop = 4];
         If[epsx >= beta1, istop = 3];
	 If[test2 <= rtol, istop = 2];
         If[test1 <= rtol, istop = 1]]];

       (* See if it is time to print something. *)
       Block[{prnt = False, debug, ww, vv, trueArnorm}, 
        If[n <= 40, prnt = True];
	If[itn <= 10, prnt = True];
        If[itn >= itnlim - 10, prnt = True];
        If[Mod[itn, 10] == 0, prnt = True];
        If[qrnorm <= 10*epsx, prnt = True];
        If[qrnorm <= 10*epsr, prnt = True];
        If[Acond <= 10^-2/eps, prnt = True];
        If[istop != 0, prnt = True];
        If[show > 1 && prnt,
	 If[Mod[itn, 10] == 0, Print[""]];
	   Print[{itn, x[[1]], test1, test2, Anorm, Acond, gbar/Anorm}]];

	debug = False (* True *);
	If[debug, (* Print true Arnorm.
                     This works only if no preconditioning. *)
         vv = b - A[x] + shift*x; (* vv=b-(A-shift*I)*x *)
         ww = A[vv] - shift*vv; (*  ww=(A-shift*I)*vv="Ar" *)
         trueArnorm = Norm[ww];
         Print["Arnorm=", Arnorm, " True||Ar|| =", trueArnorm]]];
 
       If[istop != 0, Break[]]]]]]];
  resvec = If[itn>0,Take[resvec, {itn}],0];

  (* Display final status. *)
  If[show > 0,
   Print[last, "istop=", istop, " itn=", itn];
   Print[last, "Anorm=", Anorm, " Acond=", Acond];
   Print[last, "rnorm=", rnorm, " ynorm=", ynorm];
   Print[last, "Arnorm=", Arnorm];
   Print[last, "time (seconds): M=", tm, ", A=", ta, ", orthogonalization=",
	 tortho, ", total=", SessionTime[] - tinit];
   Print[last, msg[[istop + 2]]]];
  {x, istop, itn, rnorm, Arnorm, Anorm, Acond, ynorm resvec}];
(*  End function minres.m *)

makrOrtho::usage = "Create instance of the orthogonalizer.";
makeOrtho[n_, localSizeIn_] := 
 Module[{localPointer, localOrtho, localSize = localSizeIn, 
   localVQueueFull, localV, localVEnqueue, localVOrtho},
   (* boolean to tell whether localReOrtho is on based on the value of localSize *)
  localOrtho = False;
  If[localSize > 0,
   localPointer = 0; (* tells number of prior Lanczos vectors stored *)
   localOrtho = True; (* turn on local reorthogonalization *)
   localVQueueFull = False; (* boolean that tells whether we have stored all prior localSize lanczos vectors *)
   (* Allocate storage for the number of the latest v_k^T vectors. *)
   (* Can't store more than min dimension of the matrix, n *)
   localV = Array[Null &, {n, Min[localSize, n]}]];

  (* This function stores v into the circular buffer localV *)
  localVEnqueue = 
   Function[v, 
    If[localPointer < localSize, (* localPointer counts the number currently stored *)
     localPointer = localPointer + 1, (* not full yet *)
     localPointer = 1; (* Remain orthogonal to previous localSize, so erase first one and continue circularly. *);
     localVQueueFull = True]; (* Set boolean to true for being full *)
    localV[[localPointer]] = v]; (* Store v in the column corresponding to localPointer *)

  (*Perform local reorthogonalization of v *)
  localVOrtho = 
   Function[v, 
    Block[{vOutput = v, localOrthoLimit}, 
     If[localVQueueFull, 
      localOrthoLimit = localSize, (* calculate where to terminate loop *)
      localOrthoLimit = localPointer]; (* localPointer<localSize *)
     Do[ Block[{vtemp = localV[[localOrthoCount]]},
       (* reorthogonalize 1 by 1 *)
       vOutput = vOutput - (vOutput.vtemp)*vtemp
       (* orthogonalize to each stored vector-
       note we don't have to normalize since it is explicitly done
       in the code and so we don't need to redo it *)], 
       {localOrthoCount, localOrthoLimit}];
     vOutput]];
  {localOrtho, localVEnqueue, localVOrtho}];
