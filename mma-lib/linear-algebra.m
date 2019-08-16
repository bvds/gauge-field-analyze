condition::usage = "Calculate 2-norm condition number of a matrix using the SVD."; 
condition[m_] := 
 If[Length[m] < 1000, (First[#]/Last[#]) &[
   SingularValueList[Normal[m], Tolerance -> 0]], 
  First[SingularValueList[m, 1, Method -> "Arnoldi"]/ 
	SingularValueList[m, -1, Method -> {"Arnoldi"}]]];

ersatzLanczos::usage = 
  "Perform Lanczos tri-diagonalization for a given number of steps,
then use Mathematica routine Eigensystem[] to find the eigenvalues
and eigenvectors of the tridiagonal matrix.
This is a Work-around for the fact that the Mathematica function
Eigensystem[] does not allow the Matrix to be passed as a function.
Since the resulting eigenvectors have to be transformed
back, all the Lanczos vectors must be stored.";
Options[ersatzLanczos] := {maxIterations -> Automatic, 
  printDetails -> 0, initialVector -> Automatic, 
  eigenPairs -> All};
(* This code is adapted from MINRES, file minres.m *)
ersatzLanczos[A_?MatrixQ, rest___] := 
    ersatzLanczos[(A.#) &, Length[A], rest]; 
ersatzLanczos[A_Function, n_Integer, M_?MatrixQ, rest___] := 
    ersatzLanczos[A, n, LinearSolve[M], rest]; 
ersatzLanczos[A_Function, n_Integer, 
	    M : (None | _Function | _LinearSolveFunction) : None, 
	    shift : _?NumberQ : 0.0, OptionsPattern[]] := 
  Block[{debug = False,
   show = Replace[OptionValue[printDetails], {False->0, True->1}],
   itnlim = 
    If[OptionValue[maxIterations] === Automatic, 4 n, 
       OptionValue[maxIterations]],
   last = "Exit ersatzLanczos.  ", 
   (* Maintain the message numbering used in MINRES.
      Some of these don't make sense for an eigensystem, *)
   msg = {"beta2=0. If M=I, b and x are eigenvectors", 
     "beta1=0.", (* 0 *) 
     "A solution to Ax=b was found, given rtol", 
     "A least-squares solution was found, given rtol" , 
     "Reasonable accuracy achieved, given machine epsilon", (* 3 *) 
     "x has converged to an eigenvector" , 
     "acond has exceeded 0.1/machine epsilon", 
     "The iteration limit was reached", (* 6 *)
     "A does not define a symmetric matrix", 
     "M does not define a symmetric matrix", 
     "M does not define a pos-def preconditioner"}, (* 9 *)
   tinit = SessionTime[], ta = 0, tm = 0, tortho = 0,
   addTime = Function[{timer, expr},
     Block[{t1 = SessionTime[], result = expr}, 
	   timer += SessionTime[] - t1; result], {HoldAll, SequenceHold}],
   istop = 0, itn = 0, done = False, beta1, vv = {},
   alpha = Array[Null&, {n}], beta = Array[Null&, {n + 1}],
   (* Set up y and v for the first Lanczos vector 
   v1.y=beta1 P^T v1,
   where P=C**(-1).v is really P^T v1. *)  
   y = If[OptionValue[initialVector] === Automatic, RandomReal[1, n],
	  OptionValue[initialVector]], r1, r2, s, v, Anorm = 0.0},

  r1 = y; (* initial guess x=0 initial residual *)
  
  If[M =!= None, addTime[tm, y = M[y]]];
  beta1 = Conjugate[r1].y;
  
  (* Test for an indefinite preconditioner.
  If b=0 exactly, stop with x=0. *)
  If[beta1 < 0, istop = 9; show = 1; done = True];
  If[beta1 == 0, show = 1; done = True]; 
  If[beta1 > 0, beta1 = Sqrt[beta1]]; (* Normalize y to get v1 later. *)

  beta[[1]] = beta1;
  r2 = r1;
  While[itn < itnlim && istop == 0, (* max num of iter *)
   
   itn = itn + 1;
   (* Obtain quantities for the next Lanczos vector vk+1,k=1,2,... 
   The general iteration is similar to the case k=1 with v0=0:
   p1=Operator*v1-beta1*v0,
       alpha1=v1^Tp1,q2=p1-alpha1*v1,
       beta2^2=q2^Tq2,
       v2=(1/beta2) q2.
   Again, y=betak P vk,where P=C**(-1). *)
   s = 1.0/beta[[itn]]; (* Normalize previous vector (in y). *)
   v = s*y; (* v=vk if P=I *)
 
   addTime[tortho, AppendTo[vv, v]];
   addTime[ta, y = A[v]];
   If[shift != 0, y -= shift*v]; (* shift is 0 otherwise solving A-shift*I *)
   If[itn >= 2, 
    y = y - (beta[[itn]]/beta[[itn - 1]])*r1]; (* normalization is the division r1 by oldb *)
   alpha[[itn]] = v.y; (* alphak *)
   y = (-alpha[[itn]]/beta[[itn]])*r2 + y; (* normalization of r2/beta=v *)

   (* Orthogonalize versus the previous lanczos vectors*)
   addTime[tortho, y = y - (vv.y).vv];

   r1 = r2; (* r1 is unnormalized vold *)
   r2 = y; (* r2 is unnormalized v *)
   If[M =!= None, addTime[tm, y = M[r2]]];
   beta[[itn + 1]] = r2.y; (* beta=betak+1^2 *)
   If[beta[[itn + 1]] < 0, istop = 9; Break[]];
   beta[[itn + 1]] = Sqrt[beta[[itn + 1]]];
   
   (* Estimate various norms. *)
   (* This estimate of Anorm was borrowed from minres1.m *)
   Block[{pnorm = 
      Norm[{alpha[[itn]], beta[[itn + 1]], 
        If[itn >= 2, beta[[itn]], 0]}], epsx, eps = $MachineEpsilon},
    Anorm = Max[Anorm, pnorm];
    epsx = Anorm*eps;
    
    (* See if any of the stopping criteria are satisfied.
       In rare cases, istop is already-1 from above (Abar=const*I). *)
    If[istop == 0,
     If[itn >= itnlim, istop = 6];
     If[epsx >= beta[[itn + 1]], istop = 3]]]
  ];
  If[debug, Print["alpha:  ", alpha];
	    Print["beta:  ", beta];
	    Print["vv: ", MatrixForm[vv]];
	    Print["Norm vv:  ", Map[Norm, vv]]];
  
  (* Display final status. *)
  If[show > 0,
   Print[last, "istop=", istop, " itn=", itn];
   Print[last, "Anorm=", Anorm];
   Print[last, "time (seconds): M=", tm, ", A=", ta, ", orthogonalization=",
	 tortho, ", total=", SessionTime[] - tinit];
   Print[last, msg[[istop + 2]]]];
  
  Block[{
      tt = SparseArray[{
	  Band[{1, 1}] -> Take[alpha, itn], 
	  Band[{2, 1}] -> Take[beta, {2, itn}], 
	  Band[{1, 2}] -> Take[beta, {2, itn}]}, {itn, itn}],
      eigenPairs = 
      If[OptionValue[eigenPairs] === All, itn, OptionValue[eigenPairs]], 
      vals, vecs}, 
	If[debug,
	   Print["T:  ", MatrixForm[tt]]; 
	   Print["A:  ", MatrixForm[Transpose[vv].tt.vv]]];
	(* Switch between dense and Arnoldi, based
	   on the number of eigenpairs sought.  *)
	(* In the Arnoldi case, we specify a starting vector
           consistent with the above tridiagonalization. *)
	{vals, vecs} = If[2 Length[Take[alpha, eigenPairs]] > itn,
	   Eigensystem[tt, eigenPairs],
 	   Eigensystem[tt, eigenPairs, Method -> {"Arnoldi",
		"StartingVector" -> Table[If[i==1,1,0],{i,itn}]}]];
	{vals, vecs.vv}]
 ];

(*<<"minres.m"

multiLinearSolve::usage = 
  "Returns a function that solves a linear system for a given
matrix A.  The method improves in speed with multiple solutions.";
Options[multiLinearSolve] := {Tolerance -> 10^-7};
multiLinearSolve[A_?MatrixQ, rest___] := 
    multiLinearSolve[(A.#) &, Length[A], rest]; 
multiLinearSolve[A_Function, _Integer, OptionsPattern[]] := 
  Module[{xx = {}, bb={}},
    Function[b, Block[{b1, x1, x2, y = bb.b, b1Norm},
      (* First, project out already-calculated directions. *)      
      b1 = b - y.bb;
      x2 = y.xx;
      b1Norm = Norm[b1];
      If[b1Norm > OptionValue[Tolerance],
	 (* Use MINRES because the Mathematica function LinearSystem[]
	    cannot handle the matrix input as a function. *)
	 x1 = minres[A, b1];
	 (* Note that bb is orthonormal *)
	 AppendTo[xx, x1/b1Norm];
	 AppendTo[bb, b1/b1Norm];
	 x1 + x2,
	 x2
      ]]]];
 *)
