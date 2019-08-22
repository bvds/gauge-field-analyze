condition::usage = "Calculate 2-norm condition number of a matrix using the SVD.";
condition[m_] :=
 If[Length[m] < 1000, (First[#]/Last[#])&[
   SingularValueList[Normal[m], Tolerance -> 0]],
  First[SingularValueList[m, 1, Method -> "Arnoldi"]/
	SingularValueList[m, -1, Method -> "Arnoldi"]]];

ersatzLanczos::usage =
  "Perform Lanczos tri-diagonalization for a given number of steps,
then use Mathematica routine Eigensystem[] to find the eigenvalues
and eigenvectors of the tridiagonal matrix.
This is a Work-around for the fact that the Mathematica function
Eigensystem[] does not allow the Matrix to be passed as a function.
Since the resulting eigenvectors have to be transformed
back, all the Lanczos vectors must be stored.
The default is to restart Lanczos after each Ritz pair is found.
Thus, less memory is needed and re-orthogonalization is speeded
up, but certainly not as fast as an implicit restart
or thick restart strategy. 
";
Options[ersatzLanczos] := {maxIterations -> Automatic,
  printDetails -> 0, printsPerDecade -> 10, initialVector -> Automatic,
  innerInterval -> 1, rTolerance -> $MachineEpsilon,
  eigenPairs -> 1, orthoSubspace -> None,
  (* These flags are associated with an attempt to restart Lanczos. *)	   
  restartLanczos -> True, anyPairs -> False};
(* This code is adapted from MINRES, file minres.m *)
ersatzLanczos[A_?MatrixQ, rest___] :=
    ersatzLanczos[(A.#)&, Length[A], rest];
ersatzLanczos[A_Function, n_Integer, M_?MatrixQ, rest___] :=
    ersatzLanczos[A, n, LinearSolve[M], rest];
ersatzLanczos[A_Function, n_Integer,
	    M:(None | _Function | _LinearSolveFunction):None,
	    shift:_?NumberQ:0.0, OptionsPattern[]] :=
  Block[{debug = False,
   show = Replace[OptionValue[printDetails], {False->0, True->1}],
   (* Since we perform full reorthogonalization, limit to n steps. *)
   itnlim = Min[n, Replace[OptionValue[maxIterations], Automatic -> n]],
   (* Choose integers so that they are evenly spaced on a log scale,
     making sure we hit powers of 10. *)
   valuesPerDecade = Floor[10^(Floor[0.5 + Log[10, #1] #2]/#2) + 0.5]&,
   first = "Enter ersatzLanczos. ",
   middlePrompt = "esatzLanczos ", 
   last = "Exit ersatzLanczos.  ",
   (* Maintain the message numbering used in MINRES.
      Some of these don't make sense for an eigensystem, *)
   msg = {"beta2=0. If M=I, b and x are eigenvectors",
     "beta1=0.", (* 0 *)
     "Requested eigenpairs were found, given rtol",
     "A least-squares solution was found, given rtol",
     "Reasonable accuracy achieved, given machine epsilon", (* 3 *)
     "x has converged to an eigenvector",
     "acond has exceeded 0.1/machine epsilon",
     "The iteration limit was reached", (* 6 *)
     "A does not define a symmetric matrix",
     "M does not define a symmetric matrix",
     "M does not define a pos-def preconditioner"}, (* 9 *)
   tinit = SessionTime[], ta = 0, tm = 0, tortho = 0, tinner = 0,
   addTime = Function[{timer, expr},
     Block[{t1 = SessionTime[], result = expr},
	   timer += SessionTime[] - t1; result], {HoldAll, SequenceHold}],
   istop = 0, itn = 0, done = False, vv = {}, alpha, beta,
   (* Set up y and v for the first Lanczos vector
   v1.y=beta1 P^T v1,
   where P=C**(-1).v is really P^T v1. *)
   y = If[OptionValue[initialVector] === Automatic,
	  RandomReal[{-1.0, 1.0}, n], OptionValue[initialVector]],
   r1 = Null, r2, v, Anorm = 0.0,
   rnorm = Null, vals = Null, vecs = Null, beigen,
   take = OptionValue[eigenPairs]},
  alpha = Array[Null&, {itnlim}];
  beta = Array[Null&, {itnlim + 1}];
  (* Orthogonalize initial vector with respect to some subspace. *)
  If[OptionValue[orthoSubspace] =!= None,
     y = OptionValue[orthoSubspace][y]];
  r2 = y; (* initial guess x=0 initial residual *)

  (* Procedure that solves the inner eigenvalue problem,
    setting values for vals, vecs, and rnorm.. *)
  beigen = Function[Block[{
      tt = SparseArray[{
	  Band[{1, 1}] -> Take[alpha, #],
	  Band[{2, 1}] -> Take[beta, {2, #}],
	  Band[{1, 2}] -> Take[beta, {2, #}]}, {#, #}]},
	If[debug,
	   Print["T:  ", MatrixForm[tt]];
	   Print["A:  ", MatrixForm[Transpose[vv].tt.vv]]];
	{vals, vecs} = Eigensystem[tt, take, Method -> "Banded"]];
	(* Convergence test from
          http://www.netlib.org/utk/people/JackDongarra/etemplates/node103.html *)
	rnorm = Abs[beta[[# + 1]] Map[Last, vecs]]];

  If[show>0,
     Print[first,"n=",n," itnlim=", itnlim,
	   " rTolerance=", OptionValue[rTolerance]]];

  If[M =!= None, addTime[tm, y = M[y]]];
  beta[[1]] = Conjugate[r2].y;

  (* Test for an indefinite preconditioner.
  If b=0 exactly, stop with x=0. *)
  If[beta[[1]] < 0, istop = 9; show = 1; done = True];
  If[beta[[1]] == 0, show = 1; done = True];
  If[beta[[1]] > 0,
     beta[[1]] = Sqrt[beta[[1]]]]; (* Normalize y to get v1 later. *)

  While[itn < itnlim && istop == 0, (* max num of iter *)
 
   itn = itn + 1;
   (* Obtain quantities for the next Lanczos vector vk+1,k=1,2,...
   The general iteration is similar to the case k=1 with v0=0:
   p1=Operator*v1-beta1*v0,
       alpha1=v1^Tp1,q2=p1-alpha1*v1,
       beta2^2=q2^Tq2,
       v2=(1/beta2) q2.
   Again, y=betak P vk,where P=C**(-1). *)
   Block[{s = 1.0/beta[[itn]]}, (* Normalize previous vector (in y). *)
	 v = s*y]; (* v=vk if P=I *)

   addTime[tortho, AppendTo[vv, v]];
   addTime[ta, y = A[v]];
   If[shift != 0, y -= shift*v]; (* shift is 0 otherwise solving A-shift*I *)
   If[itn >= 2,
      y = y - (beta[[itn]]/beta[[itn - 1]])*r1]; (* normalization is the division r1 by oldb *)
   alpha[[itn]] = v.y; (* alphak *)
   y = y - (alpha[[itn]]/beta[[itn]])*r2; (* normalization of r2/beta=v *)

   (* Orthogonalize versus the previous lanczos vectors
    using modified Gram-Schmidt.
    Using Scan[...] is much slower than Do[...]. *)
   addTime[tortho, Do[y = y - (vv[[i]].y) vv[[i]], {i, Length[vv]}]];

   (* Orthogonalize with respect to some subspace. *)
   If[OptionValue[orthoSubspace] =!= None,
      y = OptionValue[orthoSubspace][y]];
 
   r1 = r2; (* r1 is unnormalized vold *)
   r2 = y; (* r2 is unnormalized v *)
   If[M =!= None, addTime[tm, y = M[r2]]];
   beta[[itn + 1]] = r2.y; (* beta=betak+1^2 *)
   If[beta[[itn + 1]] < 0, istop = 9; Break[]];
   beta[[itn + 1]] = Sqrt[beta[[itn + 1]]];

   (* Periodically solve inner eigenvalue problem so we can
      estimate convergence. *)
   If[itn > 1 && itn >= Abs[take] &&
      Mod[itn, OptionValue[innerInterval]] == 0,
      addTime[tinner, beigen[itn]];
      If[Max[rnorm] <= OptionValue[rTolerance],
	 istop = 1,
	 (* This case is associated with the restart strategy *)
	 If[OptionValue[anyPairs] &&
	    Min[rnorm] <= OptionValue[rTolerance],
	    istop = 4]],
      {vals, vecs} = {Null, Null}];

   (* Estimate various norms. *)
   (* This estimate of Anorm was borrowed from minres1.m *)
   Block[{pnorm =
      Norm[{alpha[[itn]], beta[[itn + 1]],
        If[itn >= 2, beta[[itn]], 0]}], epsx, eps = $MachineEpsilon},
    Anorm = Max[Anorm, pnorm];
    epsx = Anorm*eps;
  
    (* See if any of the stopping criteria are satisfied.
       In rare cases, istop is already -1 from above (Abar=const*I). *)
    If[istop == 0,
     If[itn >= itnlim, istop = 6];
     If[epsx >= beta[[itn + 1]], istop = 3]]];

   If[show > 1 && itn == valuesPerDecade[itn, OptionValue[printsPerDecade]],
      Print[middlePrompt, " itn=", itn,
	    " rnorm=", {Min[rnorm], Max[rnorm]},
	    " time=", SessionTime[] - tinit]]

  ];  (* end of main loop *)

  If[debug, Print["alpha:  ", alpha];
	    Print["beta:  ", beta];
	    Print["vv: ", MatrixForm[vv]];
	    Print["Norm vv:  ", Map[Norm, vv]]];

  If[vals === Null, beigen[itn]];

  (* Filter out converged eigenpairs.
    This will only have an effect when istop != 1. *)
  Block[{picker = Map[(# < OptionValue[rTolerance])&, rnorm]},
	vals = Pick[vals, picker];
	vecs = Pick[vecs, picker]];

  (* Display final status. *)
  If[show > 0,
     Print[last, "istop=", istop, " itn=", itn,
	   "  converged eigenpairs=", Length[vals]];
     Print[last, "Anorm=", Anorm, " rnorm=", {Min[rnorm], Max[rnorm]}];
     Print[last, "time (seconds): M=", tm, ", A=", ta,
	   ", orthogonalization=", tortho, " inner system=", tinner,
	   ", total=", SessionTime[] - tinit];
     Print[last, msg[[istop + 2]]]];

  {vals, If[Length[vecs] > 0, vecs.vv, vecs]}
  ]/; OptionValue[anyPairs] || Not[OptionValue[restartLanczos]];

(* This is an attempt at a quick and dirty restart scheme.
  However, it is considerably slower than the original calculation. *)
ersatzLanczos[A_Function, n_Integer,
	    M:(None | _Function | _LinearSolveFunction):None,
	    shift:_?NumberQ:0.0, opts:OptionsPattern[]] :=
Block[{project, i=0,  va, ve},
Module[{vals = {}, vecs = {}},
    While[
	Length[vals] < Abs[OptionValue[eigenPairs]] &&
	i < Abs[OptionValue[eigenPairs]],
	i += 1;
	project = Function[x,
	    Block[{y = If[OptionValue[orthoSubspace] =!= None,
			  OptionValue[orthoSubspace][x], x]},
		  (* Perform modified Gram-Schmidt. *)
		  Do[y = y - (vecs[[i]].y) vecs[[i]],
		     {i, Length[vecs]}]; y]];

	{va, ve} = ersatzLanczos[
	    A, n, M, shift,
	    anyPairs -> True,
	    eigenPairs -> Sign[OptionValue[eigenPairs]] *
			(Abs[OptionValue[eigenPairs]] - Length[vals]),
	    orthoSubspace -> project,
	    opts];
	If[Length[va] == 0, Break[]];
	vals = Join[vals, va];
	vecs = Join[vecs, ve]
    ];
    {vals, vecs}
]]/; Not[OptionValue[anyPairs]] && OptionValue[restartLanczos];
