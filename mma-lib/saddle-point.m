2(*
  Define a basis of lattice degrees of freedom.

  Define operators on that basis for infinitesimal
  gauge transforms, the gradient, and the Hessian.

  Perform a search step for finding the nearest saddle point.
*)


(* Random utility functions. *)

methodName[x_] := If[Head[x] === List, First[x], x];
(* Tricky:  explicit Sequence[...] will disappear! *)
methodOptions[x_] := Apply[Sequence,
			   If[Head[x] === List, Rest[x], {}]];
SetAttributes[addTime, {HoldAll, SequenceHold}];
addTime[timer_, expr_] :=
 Block[{t1 = SessionTime[], result = expr},
  timer += SessionTime[] - t1; result];
printMemory[n_] := EngineeringForm[N[ByteCount[n]], 3];


(* Shift cutoff *)

applyCutoff1::usage = "Rescale shifts for a single link.  Large zzz \
removes the cutoff.  When Norm[shift]=Pi, there is an inflection \
point in the associated color matrix meaning that the quadratic \
approximation becomes qualitatively different than the full matrix \
exponential.  Thus, we set shift to zero when Norm[shift]>=Pi.  The \
cutoff represents the Norm[shift] where the quadratic approximation \
*starts* breaking down.  In that case, we scale back the size of the \
shift.";
applyCutoff1[hess_?VectorQ, grad_, cutoff_: 1, zzz_: 1] :=
 Block[ (* First, see if any component exceeds the cutoff.  Thus,
  we avoid any divide-by-zero error. *)
     {tooBig = Inner[(Pi zzz Abs[#1] <= Abs[#2])&, hess, grad, Or], shift},
  If[tooBig, Map[0&, grad],
   (* Otherwise, rescale based on the norm of the shift *)
   shift = grad/hess;
   Which[Pi zzz <= Norm[shift], Map[0&, shift],
    Norm[shift] < cutoff  zzz, shift,
    True, (cutoff/Norm[shift]) shift ]]];
shiftNorm[shifts_] := Max[Map[Norm, Partition[shifts, nc^2 - 1]]];
compareVectors[a_, b_] := {{shiftNorm[a], shiftNorm[b]}, {Norm[a], Norm[b]},
  a.b/(Norm[a] Norm[b])};
applyCutoff2::usage = "Rescale shifts on an entire lattice so that \
the largest norm of a shift on a link is less than the cutoff.";
applyCutoff2[shifts_?ArrayQ, cutoff_: 1, zzz_: 1] :=
 Block[{maxNorm = shiftNorm[shifts]},
  If[maxNorm < cutoff zzz, shifts,
   Print["applyCutoff2: rescale maxNorm ", maxNorm, " to ",
    cutoff]; (cutoff/maxNorm) shifts]];
applyCutoff3::usage = "Apply cutoff to the eigenspace of the Hessian. \
 In this case, the shifts are independent, but have no direct \
meaning. We infer the effect on the lattice links by using the \
rotation back onto the lattice.  We demand that the norm of the \
largest shift on any single link be less than the cutoff (default \
value Pi).";
applyCutoff3[hess_, grad_, proj_, cutoff_: Pi, zzz_: 1] :=
 Block[{maxLinkShift = Map[shiftNorm, proj], result,
	count = 0, firstValue = Null, lastValue = Null, doPairs = False},
  If[doPairs, goodPairs0 = {}; badPairs0 = {}];
  result = MapThread[
      If[cutoff zzz Abs[#1] <= Abs[#2] #3,
	 count += 1;
	 If[firstValue === Null, firstValue = #4];
	 lastValue = #4;
	 If[doPairs, AppendTo[badPairs0, {#1, #2}]]; 0,
	 If[doPairs, AppendTo[goodPairs0, {#1, #2}]]; #2/#1]&,
	{hess, grad, maxLinkShift, Range[Length[grad]]}];
  If[count > 0,
     Print["applyCutoff3:  ", count, " of ", Length[hess],
	   " zeros, between ", {firstValue, lastValue}]];
  result];
cutoffNullspace::usage =
  "Select vectors in proj that are associated with shifts that \
violate the cutoff.";
cutoffNullspace[hess_, grad_, proj_, cutoff_:Pi, zzz_:1] :=
  Block[{result, count = 0, firstValue = Null, lastValue = Null},
   result = MapThread[
       If[zzz < Infinity && cutoff zzz Abs[#1] <= Abs[#2] shiftNorm[#3],
          count += 1;
	  If[firstValue === Null, firstValue = #4];
	  lastValue = #4;
	  #3,
	  Nothing]&,
	 {hess, grad, proj, Range[Length[grad]]}];
   If[count > 0 || True,
      Print["cutoffNullSpace:  ", count, " of ", Length[hess],
	    " zeros, between ", {firstValue, lastValue}]];
   result];


(*
Define basis for describing lattice-wide shifts.   The parameter "fixed"
describes the type of shift:
    fixed == -1:  All possible shifts, including those associated with
                  gauge transformations.
    fixed == 0:   Shifts on a smaller basis where all shifts associated
                  with gauge transformations have been removed.
    fixed > 0:    For each Polyakov loop in the direction = fixed, apply
                  a constant shift.  This is compatible with choosing an
                  Axial gauge in that direction
*)

nGrad::usage = "Size of basis for lattice-wide shifts.";
nGrad[fixed_] := If[fixed < 1, nd,
   nd - 1 + 1/latticeDimensions[[fixed]]] latticeVolume[]*(nc^2 - 1);
coords2grad::usage = "Use latticeIndex, rather than linearSiteIndex \
to order the sites, so we can handle fixed>0.";
coords2grad[fixed_][dir_, coordsIn_, generator_] :=
 generator + (nc^2 - 1)*
   Block[{coords = coordsIn, dimensions = latticeDimensions},
    If[fixed == dir, coords[[fixed]] = 1; dimensions[[fixed]] = 1];
    latticeIndex[coords, dimensions] - 1 +
     latticeVolume[]*
      Sum[If[fixed == k, 1/latticeDimensions[[k]], 1], {k, dir - 1}]];
 nGauge::usage = "Size of basis for a gauge tranformation.  fixed=0: \
arbitrary gauge tranformation.  fixed>0:  gauge transforms that \
preserve Axial gauge in the \"fixed\" direction.";
nGauge[fixed_] := latticeVolume[] (nc^2 - 1) If[fixed < 1, 1,
    1/latticeDimensions[[fixed]]];
coords2gauge[fixed_][coordsIn_, generator_] :=
  generator + (nc^2 - 1)*
    Block[{coords = coordsIn, dimensions = latticeDimensions},
     If[fixed > 0, coords[[fixed]] = 1; dimensions[[fixed]] = 1];
     latticeIndex[coords, dimensions] - 1];
gradLinkTake[fixed_][dir_, coords_] :=
    {1, nc^2 - 1} + coords2grad[fixed][dir, coords, 0];
symAdd::index = "Invalid index.";
SetAttributes[symAdd, HoldFirst];
symAdd[m_, i_, j_, value_] := (
    m[{i, j}] = Lookup[m, Key[{i, j}], 0.0] + value;
    m[{j, i}] = Lookup[m, Key[{j, i}], 0.] + value);
SetAttributes[oneAdd, HoldFirst];
oneAdd[m_, i_, j_, value_] := m[{i, j}] = Lookup[m, Key[{i, j}], 0.0] + value;
reTrDot::usage = "Compute Re[Tr[a.b]] where a and b are color matrices.";
reTrDot = Compile[{{x, _Complex, 2}, {y, _Complex, 2}},
	Sum[Re[x[[i, j]] y[[j, i]]], {i, nc}, {j, nc}], {{nc, _Integer}}];


(* Gauge Transform Operator *)

gaugeTransformShifts::usage =
  "Return a matrix containing shifts associated with all possible \
infinitesimal gauge transformations.  The result is returned as a \
SparseArray.";
gaugeTransformShifts::axial = "Not Axial gauge.";
SetAttributes[gaugeTransformShifts, HoldFirst];
gaugeTransformShifts[rootGaugeField_, fixed_: - 1] :=
 ParallelSum[
  Block[{gaugeShifts = Association[], gen = suGenerators[],
    c2grad = coords2grad[fixed]},
   Do[Block[{coords = latticeCoordinates[i], gindex, index1, index2,
      u1, u2, ugu1, ugu2}, u1 = getLink[rootGaugeField][dir, coords];
     u2 = getLink[rootGaugeField][dir, shift[dir, coords, -1]];(*
     Only do the axial direction once.
     But verify that lattice is in Axial gauge. *)
     If[dir == fixed && coords[[fixed]] > 1,
	If[AnyTrue[Abs[u1 - u2], (# > 10^-6)&, 2],
	   Message[gaugeTransformShifts::axial]; Return[$Failed]];
	Continue[]];
     gindex = coords2gauge[fixed][coords, gac];
     index1 = c2grad[dir, coords, 0];
     ugu1 = ConjugateTranspose[u1].gen[[gac]].u1;
     index2 = c2grad[dir, shift[dir, coords, -1], 0];
     ugu2 = u2.gen[[gac]].ConjugateTranspose[u2];
     Do[oneAdd[gaugeShifts, gindex, index1 + ac,
       2 reTrDot[ugu1, gen[[ac]]]];
      oneAdd[gaugeShifts, gindex,
       index2 + ac, -2 reTrDot[ugu2, gen[[ac]]]], {ac,
       nc^2 - 1}]], {i, kernel, latticeVolume[], $KernelCount}, {gac,
     nc^2 - 1}, {dir, nd}];
   SparseArray[
    Normal[gaugeShifts], {nGauge[fixed],
     nGrad[fixed]}]], {kernel, $KernelCount}];


(* Hessian and gradient *)

Options[latticeHessian] = {allLinks -> True, fixedDir -> -1};
latticeHessian::usage = "Set allLinks to False to compare with \
single-link code.  Option fixedDir>0 will apply a constant shift to \
any link a given Polyakov loop winding in that direction, compatible \
with choosing an Axial gauge in that direction.";
latticeHessian[getRootLink_, OptionsPattern[]] :=
 Block[{fixed = OptionValue[fixedDir]},
 (* Adding elements to a SparseArray one at a time is very inefficient
    in Mathematica; see
    https://mathematica.stackexchange.com/questions/777/efficient-by-element-updates-to-sparsearrays.
    Instead, we accumulate Array elements in an Association and create a
    SparseArray at the end. *)
  Block[{hess = Association[], grad = Array[0.0&, nGrad[fixed]],
    gen = suGenerators[] + 0.0 I,
    sym = 0.5*Outer[(#1.#2 + #2.#1)&, suGenerators[], suGenerators[], 1] +
	  0.0 I,
    coords},
   innerLoopTime = 0;
   Do[coords = latticeCoordinates[i];
    (* Number of color matrix multiplications:
    (number of plaquettes) (16 + nc (16 + 10 nc)) *)
    Block[{
      c00 = ConjugateTranspose[getRootLink[dir2, coords]].getRootLink[dir1,
          coords],
      c10 = getRootLink[dir1, coords].getRootLink[dir2,
         shift[dir1, coords]],
      c11 = getRootLink[dir2, shift[dir1, coords]].ConjugateTranspose[
         getRootLink[dir1, shift[dir2, coords]]],
      c01 = ConjugateTranspose[
         getRootLink[dir1, shift[dir2, coords]]].ConjugateTranspose[
         getRootLink[dir2, coords]], s1, s2, s3, s4, v00, v10, v11,
      v01, w1, w2, w3, w4},
     s1 = c00.c10; s2 = c10.c11; s3 = c11.c01; s4 = c01.c00;
     v00 = c01.s1 (* =s4.c10 *);
     v10 = c00.s2 (* =s1.c11 *);
     v11 = c10.s3 (* =s2.c01 *);
     v01 = c11.s4 (* =s3.c00 *);
     w1 = s4.s2; w2 = s1.s3; w3 = s2.s4; w4 = s3.s1;
     Do[Block[{g1 = coords2grad[fixed][dir1, coords, ca1],
        g2 = coords2grad[fixed][dir2, shift[dir1, coords], ca1],
        g3 = coords2grad[fixed][dir1, shift[dir2, coords], ca1],
        g4 = coords2grad[fixed][dir2, coords, ca1],
        ss42 = s4.gen[[ca1]].s2, ss13 = s1.gen[[ca1]].s3,
        vv0011 = v00.gen[[ca1]].c11, vv1001 = v10.gen[[ca1]].c01,
        vv1100 = v11.gen[[ca1]].c00, vv0110 = v01.gen[[ca1]].c10, t0, z},
       (* Negative for inverse links. *)
       grad[[g1]] += -Im[Tr[w3.gen[[ca1]]]];
       grad[[g3]] += Im[Tr[w1.gen[[ca1]]]];
       grad[[g2]] += -Im[Tr[w4.gen[[ca1]]]];
       grad[[g4]] += Im[Tr[w2.gen[[ca1]]]];
       (* Time is dominated by matrix construction,
       and matrix construction is dominated by the inner color loop. *)
       t0 = SessionTime[];
       Do[Block[{gg1 = g1 - ca1 + ca2, gg2 = g2 - ca1 + ca2,
          gg3 = g3 - ca1 + ca2, gg4 = g4 - ca1 + ca2},
         oneAdd[hess, g1, gg1, -reTrDot[sym[[ca1, ca2]], w3]];
         oneAdd[hess, g3, gg3, -reTrDot[sym[[ca1, ca2]], w1]];
         oneAdd[hess, g2, gg2, -reTrDot[sym[[ca1, ca2]], w4]];
         oneAdd[hess, g4, gg4, -reTrDot[sym[[ca1, ca2]], w2]];
         If[OptionValue[allLinks],
          symAdd[hess, g1, gg3, reTrDot[ss42, gen[[ca2]]]];
          symAdd[hess, g2, gg4, reTrDot[ss13, gen[[ca2]]]];
          symAdd[hess, g2, gg3, reTrDot[vv0011, gen[[ca2]]]];
          symAdd[hess, g3, gg4, -reTrDot[vv1001, gen[[ca2]]]];
          symAdd[hess, g4, gg1, reTrDot[vv1100, gen[[ca2]]]];
          symAdd[hess, g1, gg2, -reTrDot[vv0110, gen[[ca2]]]]]],
	  {ca2, nc^2 - 1}]; innerLoopTime += SessionTime[] - t0],
	{ca1, nc^2 - 1}]],
      {dir1, nd - 1}, {dir2, dir1 + 1, nd}, {i, latticeVolume[]}];
   {SparseArray[Normal[hess]], grad}]];

minresLabels = Association["MINRES" -> minres, "MINRES1" -> minres1,
  "MINRES-QLP" -> minresqlp0, "MINRES-QLP0" -> minresqlp0];

dynamicPart::usage = "Dynamic (non gauge shift) part of a shift vector.
This returns a function that can be applied to a shift vector.";
dynamicPart::method = "Unknown method.";
Options[dynamicPart] = {Method -> "MINRES", printDetails -> True};
dynamicPart[gaugeTransformShifts_, OptionsPattern[]] :=
 Block[{
    printer = If[OptionValue[printDetails],
      Print["dynamicPart memory:  ", printMemory[#],
       ", gaugeTransform shifts ", printMemory[gaugeTransformShifts],
       " (bytes) for ", methodName[OptionValue[Method]]]]&},
  Which[
   OptionValue[Method] === "Orthogonalize",
   Module[{ortho = SparseArray[Chop[Orthogonalize[gaugeTransformShifts]]]},
    printer[ortho]; (# - (ortho.#).ortho)&],
   methodName[OptionValue[Method]] === "LeastSquares",
   (* Since the gauge fields were calculated in single precision,
   we choose a compatible default tolerance. 
   This tends to be slow, since the matrix is poorly conditioned.
   Maybe use some kind of restart strategy? *)
   (* Default to sparse method, using Krylov, possibly QMR? *)
   Module[{projection =
      LeastSquares[Transpose[gaugeTransformShifts], #,
        methodOptions[OptionValue[Method]], Method -> "Krylov",
        Tolerance -> 10.0^-7]&},
    printer[projection]; (# - projection[#].gaugeTransformShifts)&],
   (* In the dense case, this will be the fastest since the
      matrix factorization is pre-computed. 
      If Krylov and "Method" -> "ConjugateGradient" are selected,
      this is quite slow for some unknown reason. *)
   methodName[OptionValue[Method]] === "LinearSolve",
   Module[{linearSolve =
      LinearSolve[
       gaugeTransformShifts.Transpose[gaugeTransformShifts],
       methodOptions[OptionValue[Method]]]},
    printer[linearSolve]; (# -
       linearSolve[gaugeTransformShifts.#].gaugeTransformShifts)&],
   KeyExistsQ[minresLabels, methodName[OptionValue[Method]]],
   (* Store options, setting any default values.  rTolerance=
     10^-13 produces < 10^-7 errors in the projection. *)
   Module[{opts = {methodOptions[OptionValue[Method]],
		   rTolerance -> 10^-13},
	   minresFun = minresLabels[methodName[OptionValue[Method]]]},
      printer[Null];
      Function[x, x - First[minresFun[
          (gaugeTransformShifts.(#.gaugeTransformShifts))&,
	  gaugeTransformShifts.x, Apply[Sequence, opts]]].
			   gaugeTransformShifts]],
   True,
   Message[dynamicPart::method]; $Failed]];


(* Find saddle point *)

(*
   Find the solution to an eigensystem subject to a
   linear constraint.

   This is the computational task discussed in:
   "Large sparse symmetric eigenvalue problems
   with homogeneous linear constraints: the
   Lanczos process with inner–outer iterations"
   Gene H. Golub, Zhenyue Zhang, Hongyuan Zha
   https://core.ac.uk/download/pdf/81957508.pdf
 *)

findDelta::usage = "Given a Hessian matrix and a gradient, find the
associated shifts using either dense or sparse matrix techniques.";
findDelta::lastPairs = "Mathematica orders eigenvalues/vectors in
order of decreasing magnitude.  Thus, for large shifts, one must find
the last eigenvalues/vectors.";
findDelta::external = "Error in external program; see `1` for details.";
Options[findDelta] = {dynamicPartMethod -> Automatic, printDetails -> False,
  Method -> Automatic, rescaleCutoff -> 1, dampingFactor -> 1,
  storePairs -> False, storeHess -> False, largeShiftOptions -> {},
  storeBB -> False, debugProj -> False, largeShiftCutoff -> 2,
  (* Roughly speaking, the relative error in the matrix norm of
     the link goes as Norm[shift]^3/48. *)
  cutoffValue -> 2};

(* Dense matrix version (default) *)
findDelta[hessIn_, gradIn_, gaugeTransformShifts_, OptionsPattern[]] :=
  Block[{hess, grad, zzz = OptionValue[rescaleCutoff],
    cutoff = OptionValue[cutoffValue],
    damping = OptionValue[dampingFactor], values, oo, proj, shifts},
    If[MatrixQ[gaugeTransformShifts],
       (* Find projection onto subspace orthogonal to all gauge transforms.
       proj is dense. *)
       proj = NullSpace[gaugeTransformShifts];
       hess = proj.hessIn.Transpose[proj];
       grad = proj.gradIn,
       (* Alternatively, no projection onto the subspace. *)
       proj = IdentityMatrix[Length[hessIn], SparseArray];
       hess = hessIn;
       grad = gradIn];
    If[False, Print["Dynamic hess.grad: ",(hess.grad).proj]];
    {values, oo} = Eigensystem[Normal[hess]];
    shifts = applyCutoff3[values, oo.grad, oo.proj,
	OptionValue[largeShiftCutoff], zzz].oo.proj;
    Print["shifts norms and rescale:  ",
	  {shiftNorm[shifts], Norm[shifts],
	   shifts.hessIn.shifts/shifts.gradIn}];
    If[OptionValue[storeHess],
       hess0 = hess; grad0 = grad; proj0 = proj; oo0 = oo];
    If[OptionValue[storePairs],
       pairs0 = Transpose[{values, oo.grad}];
       shifts0 = shifts];
    -damping applyCutoff2[shifts, cutoff, zzz]] /;
  OptionValue[Method] === Automatic || OptionValue[Method] === "Dense";

(* MINRES/MINRES-QLP version. *)
findDelta[hess_, grad_, gaugeTransformShifts_, opts:OptionsPattern[]] :=
 Block[{zzz = OptionValue[rescaleCutoff],
    minresFun = minresLabels[methodName[OptionValue[Method]]],
    cutoff = OptionValue[cutoffValue],
    damping = OptionValue[dampingFactor], result, shifts, dp, dp0,
    smallProj, smallProj0, tinit = SessionTime[],
    tdp = 0, tsp = 0, tdot = 0},
   dp0 = If[MatrixQ[gaugeTransformShifts],
	    dynamicPart[gaugeTransformShifts,
			Method -> OptionValue[dynamicPartMethod]],
	    Identity];
    If[False, Print["Dynamic hess.grad: ", dp0[hess.grad]]];
   Clear[bb0]; bb0 = {};
   dp = ((If[OptionValue[storeBB] == True, AppendTo[bb0, #]];
	 addTime[tdp, dp0[#]])&);
   (* Compare sparse strategy against dense projection. *)
   If[OptionValue[debugProj],
     Module[{projDense = Nullspace[gaugeTransformShifts]},
       projDense = Transpose[projDense].projDense;
       dp = (Block[{p0=projDense.#, p1=dp0[#]},
	 If[Norm[p0-p1]>10^-10,
	     Print["** Dp discrepency: ", Norm[p0-p1]]; Abort[]]; p1]&)]];
   (* Calculate eigenvectors associated with shifts that are too
     large.  Since, the Mma routine Eigensystem[...] doesn't allow
     the matrix to be given as a function, we use ersatzLanczos[]
     as a work-around. *)
   smallProj0 = Block[{vals, vecs,
      options = OptionValue[largeShiftOptions]},
     (* Make sure we are looking for the smallest eigenvalues *)
     If[NumberQ[eigenPairs/.options] && (eigenPairs/.options)>0,
	Message[findDelta::lastPairs]; Return[$Failed]];
     {vals, vecs} = ersatzLanczos[
	 dp[addTime[tdot, hess.#]]&, Length[grad],
	 Apply[Sequence, options],
	 (* Need to determine a reasonable estimate for this. *)
	 eigenPairs -> -Min[120, Length[grad]],
         orthoSubspace -> dp,
	 initialVector -> grad, printDetails -> 1];
     Module[{largeShiftSpace = cutoffNullspace[vals, vecs.grad, vecs,
         OptionValue[largeShiftCutoff], OptionValue[rescaleCutoff]]},
	    If[Length[largeShiftSpace] > 0,
	       (# - (largeShiftSpace.#).largeShiftSpace)&, Identity]]];
   smallProj = (addTime[tsp, smallProj0[#]]&);
   (* Call the appropriate MINRES routine, setting some default values. *)
   result =
    minresFun[smallProj[dp[addTime[tdot, hess.#]]]&,
     smallProj[dp[grad]], methodOptions[OptionValue[Method]],
     Apply[Sequence,
      FilterRules[{localSize -> 5000, printDetails -> 1},
       Options[minresFun]]]]; shifts = result[[1]];
   (* Verify result is in the subspace *)
   If[False,
      Print["Result overlap with small shift subspace, \
            dynamicPart, and both:  ",
	    {shifts.smallProj[shifts], shifts.dp[shifts],
	     shifts.smallProj[dp[shifts]]}/shifts.shifts]];
   Print["findDelta time (seconds):  dynamic part=", tdp,
    ", small shifts=", tsp, ", A=", tdot, ", total=",
    SessionTime[] - tinit];
   If[OptionValue[storeHess],
      hess0 = hess; grad0 = grad; proj0 = Null; oo0 = Null];
   If[OptionValue[storePairs],
      pairs0 = Null;
      shifts0 = shifts];
   -damping applyCutoff2[shifts, cutoff, zzz]] /;
 KeyExistsQ[minresLabels, methodName[OptionValue[Method]]];

(* External version.
   As a first attempt, dump to files and run program.
   Later, maybe try LibraryLink or MathLink/WSTP.

https://mathematica.stackexchange.com/questions/199925/can-we-link-mathematica-and-fortran-with-wstp
https://mathematica.stackexchange.com/questions/6325/mathematica-library-link-how-to-use-non-standard-mint-e-g-unsigned-int-or
https://stackoverflow.com/questions/6537457/gcc-installed-mathematica-still-wont-compile-to-c
https://mathematica.stackexchange.com/questions/31545/returning-multiple-results-from-a-librarylink-function
https://mathematica.stackexchange.com/questions/8438/minimal-effort-method-for-integrating-c-functions-into-mathematica

https://reference.wolfram.com/language/LibraryLink/ref/callback/AbortQ.html

https://github.com/DaveGamble/cJSON
 *)
findDelta[hess_, grad_, gaugeTransformShifts_, opts:OptionsPattern[]] :=
 Block[{zzz = OptionValue[rescaleCutoff],
    cutoff = OptionValue[cutoffValue],
    damping = OptionValue[dampingFactor], shifts,
    tinit = SessionTime[],
    sparseExport = Function[MapThread[
        Append,
        (* Switch to zero-based array indexing *)
        {#["NonzeroPositions"] - 1, #["NonzeroValues"]}]],
    symbolString = (a_Symbol -> b_) :> SymbolName[a] -> b,
    outFile = "hess-grad-gauge.json", out, logFile = "shifts.log"},
   (* Dump dimensions, constants, and options into JSON file.
     Dump matrices and vectors: 
         hess, grad, gaugtransformationShifts
     The methods themselves will be chosen at compile time,
     since we can always compare with the mathematica value.
    *)
   Export[outFile, {
       "nc" -> nc,
       "n" -> Length[grad],
       "rescaleCutoff" -> zzz,
       "largeShiftCutoff" ->
           OptionValue[largeShiftCutoff]/.Infinity -> 10.0^20,
       "linearSolveOptions" ->
	   {methodOptions[OptionValue[Method]]}/.symbolString,
       "dynamicPartOptions" ->
           {methodOptions[OptionValue[dynamicPartMethod]]}/.symbolString,
       "largeShiftOptions" -> OptionValue[largeShiftOptions]/.symbolString,
       "hessElements" -> Length[hess["NonzeroValues"]],
       "gaugeElements" -> Length[gaugeTransformShifts["NonzeroValues"]],
       "gaugeDimension" -> Length[gaugeTransformShifts]}];
   Export["hess.dat", sparseExport[hess]];
   Export["grad.dat", grad];
   Export["gauge.dat", sparseExport[gaugeTransformShifts]];
   (* Run external program *)
   out = Run["saddle-lib/shifts",
             "hess-grad-gauge.json hess.dat grad.dat gauge.dat shifts.dat",
             ">", logFile];
   If[out != 0, Message[findDelta::external, logFile]; Return[$Failed]];
   (* Read shifts from external file *)
   shifts = ReadList["shifts.dat", Number];
   Print["findDelta time (seconds):  total=", SessionTime[] - tinit];
   If[OptionValue[storeHess],
      hess0 = Null; grad0 = Null; proj0 = Null; oo0 = Null];
   If[OptionValue[storePairs],
      pairs0 = Null;
      shifts0 = shifts];
   -damping applyCutoff2[shifts, cutoff, zzz]] /;
methodName[OptionValue[Method]] === "External";

SetAttributes[applyDelta, HoldFirst];
applyDelta[tgf_, delta_, fixed_] :=
  Block[{gradMap = gradLinkTake[fixed]},
   Do[
    Block[{coords = latticeCoordinates[i]},
     tgf[[dir, linearSiteIndex[coords]]] =
       getLink[tgf][dir, coords].
       MatrixExp[I Take[delta, gradMap[dir, coords]].suGenerators[]].
       getLink[tgf][dir, coords]],
    {dir, nd}, {i, latticeVolume[]}]];

latticeSaddlePointStep::usage =
  "Setting the option fixedDir=0 will use a projection onto the \
subspace orthogonal to all gauge transforms.  This should be the most \
correct method.  Can apply this to both dense and sparse systems.";
(* For passing options through to subfunctions: 
https://mathematica.stackexchange.com/questions/353/functions-with-options/ *)
Options[latticeSaddlePointStep] =
 Join[Options[findDelta], Options[latticeHessian]];
latticeSaddlePointStep[opts:OptionsPattern[]] :=
 Block[{hess, grad, delta, tmpGaugeField = makeRootLattice[],
	t0 = SessionTime[], t1, t2, t3},
  {hess, grad} = latticeHessian[getLink[tmpGaugeField],
	Apply[Sequence, FilterRules[{opts}, Options[latticeHessian]]]];
  t1 = SessionTime[];
  (* Debug prints *)
  Which[False,
   Print["hessian difference:  ", Chop[Normal[hess] - hess0]];
   Print["gradient difference:  ", Chop[grad - gradient0]],
   False,
   Print["hessian:", Normal[hess]];
   Print["gradient:  ", grad],
   False,
   Print["hessian, first link:  ",
    Normal[Take[hess, nc^2 - 1, nc^2 - 1]]];
   Print["gradient, first link:  ", Take[grad, nc^2 - 1]]];
  delta = findDelta[hess, grad,
    If[OptionValue[fixedDir] > -1,
       gaugeTransformShifts[tmpGaugeField, OptionValue[fixedDir]],
       None],
    Apply[Sequence, FilterRules[{opts}, Options[findDelta]]]];
  (* Debug prints *)
  Which[
      False, Print["delta difference:  ", Chop[delta - delta0]],
      False, Print["delta:  ", delta],
      False, Print["delta, first link:  ", Take[delta, nc^2 - 1]]];
  t2 = SessionTime[];
  applyDelta[tmpGaugeField, delta, OptionValue[fixedDir]];
  t3 = SessionTime[];
  Print["latticeSaddlePointStep times:  ",
	{innerLoopTime, t1 - t0, t2 - t1, t3 - t2}];
  tmpGaugeField];
