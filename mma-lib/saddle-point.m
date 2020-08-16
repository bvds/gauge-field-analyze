(*
  Define a basis of lattice degrees of freedom.

  Define operators on that basis for infinitesimal
  gauge transforms, the gradient, and the Hessian.

  Perform a search step for finding the nearest saddle point.
*)


(* Random utility functions. *)

Get[DirectoryName[$InputFileName] <> "matrix.m"];
methodName[x_] := If[Head[x] === List, First[x], x];
(* Tricky:  explicit Sequence[...] will disappear! *)
methodOptions[x_] := Apply[Sequence,
			   If[Head[x] === List, Rest[x], {}]];


(* Define norms for shifts *)
(* This differs from the norm defined in the paper by
  a factor of sqrt(2), but eigenCutoffMax, eigenCutoff2, 
  and eigenCutoffRescale should match the associated lambda
  parameters in the paper. *)
shiftNormMax[shifts_] :=
    Max[Map[Norm, Partition[shifts, nc^2 - 1]]];
shiftNorm2[shifts_] :=
    Norm[shifts] Sqrt[(nc^2-1)/Length[shifts]];
compareVectors[a_, b_] := {{(shiftNormMax[a]+shiftNormMax[b])/2,
                            1-shiftNormMax[a]/shiftNormMax[b]},
                           {(Norm[a]+Norm[b])/2, 1-Norm[a]/Norm[b]},
                           Block[{cos = a.b/(Norm[a] Norm[b])},
                                 {cos, Sqrt[(1+cos) (1-cos)]}]};
matrixAsymmetry[mat_] := Max[Abs[Flatten[hess-Transpose[mat]]]];

(* Shift cutoff *)

applyCutoff1::usage = "Rescale shifts for a single link.  \
A large value for zzz removes the cutoff.  \
When Norm[shift]=Pi, there is an inflection \
point in the associated color matrix meaning that the quadratic \
approximation becomes qualitatively different than the full matrix \
exponential.  Thus, we set shift to zero when Norm[shift]>=Pi.  The \
cutoff represents the Norm[shift] where the quadratic approximation \
*starts* breaking down.  In that case, we scale back the size of the \
shift.";
applyCutoff1[hess_?VectorQ, grad_, cutoff_: Pi, zzz_: 1] :=
 Block[ (* First, see if any component exceeds the cutoff.  Thus,
  we avoid any divide-by-zero error. *)
     {tooBig = Inner[(Pi zzz Abs[#1] <= Abs[#2])&, hess, grad, Or], shift},
  If[tooBig, Map[0&, grad],
   (* Otherwise, rescale based on the norm of the shift *)
   shift = grad/hess;
   Which[Pi zzz <= Norm[shift], Map[0&, shift],
         Norm[shift] < cutoff  zzz, shift,
         True, (cutoff/Norm[shift]) shift]]];
applyCutoff2::rescale = "Rescale maxNorm `1` to `2`";
applyCutoff2::usage = "Rescale shifts on an entire lattice so that \
the largest norm of a shift on a link is less than the cutoff.";
applyCutoff2[shifts_?ArrayQ, cutoff_: 1, zzz_: 1] :=
 Block[{maxNorm = shiftNormMax[shifts]},
  If[maxNorm < cutoff zzz, shifts,
     Message[applyCutoff2::rescale, NumberForm[maxNorm, 8], cutoff];
     (cutoff/maxNorm) shifts]];
applyCutoff3::usage = "Apply cutoff to the eigenspace of the Hessian. \
 In this case, the shifts are independent, but have no direct \
meaning. We infer the effect on the lattice links by using the \
rotation back onto the lattice.  We demand that the norm of the \
largest shift on any single link be less than the cutoff (default \
value Pi).";
Options[applyCutoff3] = {
    (* When finding a minimum, (as opposed to a saddle
      point, we want to remove all negative second derivative (concave)
      components. *)
    "removeConcave" -> False,
    "eigenCutoffMax" -> N[Pi], "eigenCutoff2" -> N[Pi], "zzz" -> 1};
applyCutoff3::count = "{total, max, norm, concave} = `1` of `2` eigenpairs between `3`";
applyCutoff3[hess_, grad_, proj_, OptionsPattern[]] :=
 Block[{result, count = 0, countMax = 0, count2 = 0,
        countConcave = If[OptionValue["removeConcave"], 0, Null],
        firstValue = Null, lastValue = Null},
  result = MapThread[
      Block[{
          testConcave = OptionValue["removeConcave"] && #1 <= 0,
          testMax = OptionValue["zzz"] < Infinity &&
                    OptionValue["eigenCutoffMax"] OptionValue["zzz"] Abs[#1] <=
                               Abs[#2] shiftNormMax[#3],
          test2 = OptionValue["zzz"] < Infinity &&
                  OptionValue["eigenCutoff2"] OptionValue["zzz"] Abs[#1] <=
                             Abs[#2] shiftNorm2[#3]},
            If[testConcave, countConcave++];
            If[testMax, countMax++];
            If[test2, count2++];
            If[testMax || test2 || testConcave,
	       count += 1;
	       If[firstValue === Null, firstValue = #4];
	       lastValue = #4;
               0,
	       #2/#1]]&,
      {hess, grad, proj, Range[Length[grad]]}];
  If[count > 0,
     Message[applyCutoff3::count, {count, countMax, count2, countConcave},
             Length[hess], {firstValue, lastValue}]];
  result];
cutoffNullspace::usage =
  "Select vectors in proj that are associated with shifts that \
violate the cutoff.";
Options[cutoffNullspace] = Options[applyCutoff3];
cutoffNullspace[hess_, grad_, proj_, OptionsPattern[]] :=
 Block[{result,  count = 0, countMax = 0, count2 = 0,
        countConcave = If[OptionValue["removeConcave"], 0, Null],
        firstValue = Null, lastValue = Null},
   result = MapThread[
      Block[{
          testConcave = OptionValue["removeConcave"] && #1 <= 0,
          testMax = OptionValue["zzz"] < Infinity &&
                    OptionValue["eigenCutoffMax"] OptionValue["zzz"] Abs[#1] <=
                               Abs[#2] shiftNormMax[#3],
          test2 = OptionValue["zzz"] < Infinity &&
                  OptionValue["eigenCutoff2"] OptionValue["zzz"] Abs[#1] <=
                             Abs[#2] shiftNorm2[#3]},
            If[testConcave, countConcave++];
            If[testMax, countMax++];
            If[test2, count2++];
            If[testMax || test2,
               count += 1;
	       If[firstValue === Null, firstValue = #4];
	       lastValue = #4;
	       #3,
	       Nothing]]&,
	   {hess, grad, proj, Range[Length[grad]]}];
   If[count > 0 || True,
      Print["cutoffNullspace:  {total, max, norm, concave} = ",
            {count, countMax, count2, countConcave}, " of ", Length[hess],
	    " eigenpairs between ", {firstValue, lastValue}]];
   result];


(*
Define basis for describing lattice-wide shifts.  The parameter "fixed"
describes the type of shift:
    fixed == -1:  All possible shifts, including those associated with
                  gauge transformations.
    fixed == 0:   Shifts on a smaller basis where all shifts associated
                  with gauge transformations have been removed.
    fixed > 0:    For each Polyakov loop in the direction = fixed, apply
                  a constant shift.  This is compatible with choosing an
                  Axial gauge in that direction
*)

lattice::fixedTransBC = "The combination fixed>0 and TRANS_GAUGEBC not implemented";
nGrad::usage = "Size of basis for lattice-wide shifts.";
nGrad[fixed_] := (nc^2 - 1)*If[
    fixed < 1,
    reducedLinkCount[],
    If[latticeBC != "PERIODIC_GAUGEBC",
       Message[lattice::fixedTransBC];
       Return[$Failed]];
    latticeVolume[]*(nd - 1 + 1/latticeDimensions[[fixed]])];
coords2grad::usage = "Use latticeIndex, rather than linearSiteIndex \
to order the sites, so we can handle fixed>0.
For efficient code parallelization, links that are spatially close \
should be close in the basis.";
coords2grad[fixed_][dir_, coords_, generator_] :=
    generator + (nc^2 - 1)*
    If[fixed>0,
       If[latticeBC != "PERIODIC_GAUGEBC",
          Message[lattice::fixedTransBC];
          Return[$Failed]];
       Block[{cords = coords, dimensions = latticeDimensions},
        cords[[fixed]] = 1; dimensions[[fixed]] = 1;
        (latticeIndex[cords, dimensions] - 1) *
        (1 + (nd-1) latticeDimensions[[fixed]]) +
        If[fixed<dir, 1 + latticeDimensions[[fixed]]*(dir - 2),
           latticeDimensions[[fixed]]*(dir - 1)] +
        If[dir==fixed, 1, coords[[fixed]]] - 1],
       (* non-Axial case *)
       If[boundaryLinkQ[dir, coords],
          Infinity,
          reducedLinkIndex[dir, coords] - 1]];
nGauge::usage = "Size of basis for a gauge tranformation.  fixed=0: \
arbitrary gauge tranformation.  fixed>0:  gauge transforms that \
preserve Axial gauge in the \"fixed\" direction.";
nGauge[fixed_] := (nc^2 - 1)*
    If[fixed < 1,
       reducedSiteCount[],
       If[latticeBC != "PERIODIC_GAUGEBC",
          Message[lattice::fixedTransBC];
          Return[$Failed]];
       latiiceVolume[]/latticeDimensions[[fixed]]];
coords2gauge[fixed_][coordsIn_, generator_] :=
  generator + (nc^2 - 1)*(
   If[fixed > 0,
      Block[{coords = coordsIn, dimensions = latticeDimensions},
            If[latticeBC != "PERIODIC_GAUGEBC",
               Message[lattice::fixedTransBC];
               Return[$Failed]];
            coords[[fixed]] = 1; dimensions[[fixed]] = 1;
            latticeIndex[coords, dimensions]],
      If[boundarySiteQ[coordsIn],
         Infinity,
         reducedSiteIndex[coordsIn]]] - 1);
gradLinkTake[fixed_][dir_, coords_] :=
    {1, nc^2 - 1} + coords2grad[fixed][dir, coords, 0];


(* Gauge Transform Operator *)

gaugeTransformShifts::usage =
  "Return a matrix containing shifts associated with all possible \
infinitesimal gauge transformations.  The result is returned as a \
SparseArray.";
gaugeTransformShifts::axial = "Not Axial gauge.";
SetAttributes[gaugeTransformShifts, HoldFirst];
gaugeTransformShifts[rootGaugeField_Symbol, fixed_: -1] := (
  (* In the boundary case, make sure lookup table is initialized
   before parallelization. *)
  If[latticeBC =!= "PERIODIC_GAUGEBC",
     reducedLinkIndex[1, Table[1, {nd}]]];
 ParallelSum[
  Block[{gaugeShifts = Association[], gen = SUGenerators[],
    c2grad = coords2grad[fixed]},
   Do[Block[{coords = latticeCoordinates[i], gindex, index1, index2,
             u1, u2, ugu1, ugu2},
     u1 = getLink[rootGaugeField][dir, coords];
     u2 = getLink[rootGaugeField][dir, shift[dir, coords, -1]];
     (* Only do the axial direction once.
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
     Do[
         oneAdd[gaugeShifts, gindex, index1 + ac,
                2 reTrDot[ugu1, gen[[ac]]]];
         oneAdd[gaugeShifts, gindex,
                index2 + ac, -2 reTrDot[ugu2, gen[[ac]]]],
         {ac, nc^2 - 1}]],
      {i, kernel, latticeVolume[], $KernelCount},
      {gac, nc^2 - 1}, {dir, nd}];
   If[latticeBC =!= "PERIODIC_GAUGEBC",
      gaugeShifts = KeySelect[gaugeShifts,
                              (#[[1]]<Infinity && #[[2]]<Infinity)&]];
   SparseArray[Normal[gaugeShifts],
               {nGauge[fixed], nGrad[fixed]}]],
  {kernel, $KernelCount}]);


(* Hessian and gradient *)

Options[actionHessian] = {allLinks -> True, fixedDir -> 0,
                           fullMatrix -> False};
actionHessian::usage = "Set allLinks to False to compare with \
single-link code.  Option fixedDir>0 will apply a constant shift to \
any link a given Polyakov loop winding in that direction, compatible \
with choosing an Axial gauge in that direction.
Set the option fullMatrix to True to explicitly compute the upper triangle of the matrix.";
actionHessian[getRootLink_, OptionsPattern[]] :=
 Block[{fixed = OptionValue[fixedDir]},
  (* In the boundary case, make sure lookup table is initialized
   before parallelization. *)
  If[latticeBC =!= "PERIODIC_GAUGEBC",
     reducedLinkIndex[1, Table[1, {nd}]]];
 (* Adding elements to a SparseArray one at a time is very inefficient
    in Mathematica; see
    https://mathematica.stackexchange.com/questions/777/efficient-by-element-updates-to-sparsearrays.
    Instead, we accumulate Array elements in an Association and create a
    SparseArray at the end. *)
 ParallelSum[
  Block[{hess = Association[], grad = Array[0.0&, nGrad[fixed]],
    gen = SUGenerators[] + 0.0 I,
    sym = 0.5*Outer[(#1.#2 + #2.#1)&, SUGenerators[], SUGenerators[], 1] +
	  0.0 I,
    full = OptionValue[fullMatrix],
    coords},
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
       If[g1 < Infinity,
          grad[[g1]] += -Im[Tr[w3.gen[[ca1]]]]];
       If[g3 < Infinity,
          grad[[g3]] += Im[Tr[w1.gen[[ca1]]]]];
       If[g2 < Infinity,
          grad[[g2]] += -Im[Tr[w4.gen[[ca1]]]]];
       If[g4 < Infinity,
          grad[[g4]] += Im[Tr[w2.gen[[ca1]]]]];
       (* Time is dominated by matrix construction,
       and matrix construction is dominated by the inner color loop. *)
       t0 = SessionTime[];
       Do[Block[{
           gg1 = g1 - ca1 + ca2, gg2 = g2 - ca1 + ca2,
           gg3 = g3 - ca1 + ca2, gg4 = g4 - ca1 + ca2},
         oneAdd[hess, g1, gg1, -reTrDot[sym[[ca1, ca2]], w3], full];
         oneAdd[hess, g3, gg3, -reTrDot[sym[[ca1, ca2]], w1], full];
         oneAdd[hess, g2, gg2, -reTrDot[sym[[ca1, ca2]], w4], full];
         oneAdd[hess, g4, gg4, -reTrDot[sym[[ca1, ca2]], w2], full];
         If[OptionValue[allLinks],
            symAdd[hess, g1, gg3, reTrDot[ss42, gen[[ca2]]], full];
            symAdd[hess, g2, gg4, reTrDot[ss13, gen[[ca2]]], full];
            symAdd[hess, g2, gg3, reTrDot[vv0011, gen[[ca2]]], full];
            symAdd[hess, g3, gg4, -reTrDot[vv1001, gen[[ca2]]], full];
            symAdd[hess, g4, gg1, reTrDot[vv1100, gen[[ca2]]], full];
            symAdd[hess, g1, gg2, -reTrDot[vv0110, gen[[ca2]]], full]]],
	  {ca2, nc^2 - 1}]],
	{ca1, nc^2 - 1}]],
      {dir1, nd - 1}, {dir2, dir1 + 1, nd},
      {i, kernel, latticeVolume[], $KernelCount}];
   If[latticeBC =!= "PERIODIC_GAUGEBC",
      hess = KeySelect[hess,
                       (#[[1]]<Infinity && #[[2]]<Infinity)&]];
   {SparseArray[Normal[hess],
                {nGrad[fixed], nGrad[fixed]}], grad}],
  {kernel, $KernelCount}]];


minresLabels = Association["MINRES" -> minres, "MINRES1" -> minres1,
  "MINRES-QLP" -> minresqlp0, "MINRES-QLP0" -> minresqlp0];

dynamicPart::usage = "Dynamic (non gauge shift) part of a shift vector.
This returns a function that can be applied to a shift vector.";
dynamicPart::method = "Unknown method.";
Options[dynamicPart] = {Method -> "MINRES", printDetails -> True};
dynamicPartCounter = 0;
dynamicPart[gaugeTransformShifts_, OptionsPattern[]] :=
 Block[{
   printer = If[printLevel[OptionValue[printDetails], 3] > 2,
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
          (dynamicPartCounter++; gaugeTransformShifts.(#.gaugeTransformShifts))&,
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
findDelta::external = "External program exited with error `1`.";
findDelta::dynamicPartMethod = "Only dynamicPartMethod = Automatic is supported.";
findDelta::asymmetric = "Hessian must be symmetric unless Method->\"External\". Max asymmetry `1`";
findDelta::symmetric = "Since Method->\"External\", you can set fullMatrix -> False.";
Options[findDelta] = Join[
    Options[applyCutoff3],
    {dynamicPartMethod -> Automatic, printDetails -> True,
     Method -> "Dense", ignoreCutoff -> False, dampingFactor -> 1,
     storeHess -> False, largeShiftOptions -> {},
     storeBB -> False, debugProj -> False, externalAction -> Automatic,
     remoteHost -> "samson",
     (* For "detach" and "read", polling interval in seconds. *)
     readInterval -> 60,
     (* Used by BLAS libraries and matrixVector()
     Extensive benchmarking shows 2 threads * maximum MPI processes
     is optimal. *)         
     threads -> 2,
     (* A numerical value will invoke the MPI version.
     Number of MPI processes.*)
     processes -> Automatic,
     (* Specify the way matrices are divided among processors.
     See saddle-lib/matrix.c for details. *)
     "partitions" -> Automatic,
     "gaugePartitions" -> Automatic,
     (* Currently, hard-coded for a single AMD Epyc processor. *)
     mpiFlags -> None,
     (* Not used. *)
     chunkSize -> 1,
     (* Roughly speaking, the relative error in the matrix norm of
     one link goes as Norm[shift]^3/48. *)
     rescaleCutoff -> 2}];
SetAttributes[findDelta, HoldFirst];

(* Dense matrix version (default) *)
findDelta[{hess_, grad_, gauge_}, opts:OptionsPattern[]] :=
 Block[{hessp, gradp, zzz = If[OptionValue[ignoreCutoff], Infinity, 1],
        damping = OptionValue[dampingFactor], values, oo, proj},
    (* SymmetricMatrixQ[] has trouble with small nonzero elements *)
    If[matrixAsymmetry[hess] > 10^-12,
       Message[findDelta::asymmetric, matrixAsymmetry[hess]];
       Return[$Failed]];
    If[MatrixQ[gauge],
       (* Find projection onto subspace orthogonal to all gauge transforms.
         proj is dense. *)
       proj = NullSpace[gauge];
       hessp = proj.hess.Transpose[proj];
       gradp = proj.grad,
       (* Alternatively, no projection onto the subspace. *)
       proj = IdentityMatrix[Length[hess], SparseArray];
       hessp = hess;
       gradp = grad];
    (* Debug print to compare with call to testOp(...) in shifts.c *)
    If[False, Print["Dynamic hess.grad: ", hess.grad]];
    {values, oo} = Eigensystem[Normal[hessp]];
    stepShifts = applyCutoff3[
        values, oo.gradp, oo.proj,
        "zzz" -> zzz,
        Apply[Sequence, FilterRules[{opts}, Options[applyCutoff3]]]].oo.proj;
    Print["shifts norms and rescale:  ",
	  {shiftNormMax[stepShifts], shiftNorm2[stepShifts],
	   stepShifts.hess.stepShifts/stepShifts.grad}];
    If[OptionValue[storeHess],
       Print["Define hess0, grad0, proj0, oo0, pairs0"];
       hess0 = hess; grad0 = grad; proj0 = proj; oo0 = oo;
       pairs0 = Transpose[{values, oo.gradp}]];
    -damping applyCutoff2[stepShifts, OptionValue[rescaleCutoff], zzz]] /;
  OptionValue[Method] === "Dense";

(* Version using mathematica versions of MINRES/MINRES-QLP
   and Lanczos eigensystem solvers. *)
findDelta[{hess_, grad_, gauge_}, opts:OptionsPattern[]] :=
    Block[{zzz = If[OptionValue[ignoreCutoff], Infinity, 1],
    minresFun = minresLabels[methodName[OptionValue[Method]]],
    squareHess = Not[OptionValue["removeConcave"]],
    damping = OptionValue[dampingFactor], result, dp, dp0,
    smallProj, smallProj0, tinit = SessionTime[],
        tdp = 0, tsp = 0, tdot = 0, countdp = 0, countdot = 0, countsp = 0,
        countGauge = dynamicPartCounter},
   (* SymmetricMatrixQ[] has trouble with small nonzero elements *)
   If[matrixAsymmetry[hess] > 10^-12,
      Message[findDelta::asymmetric, matrixAsymmetry[hess]];
      Return[$Failed]];
   dp0 = If[MatrixQ[gauge],
	    dynamicPart[gauge,
			Method -> OptionValue[dynamicPartMethod]],
	    Identity];
   If[False, Print["Dynamic hess.grad: ", dp0[hess.grad]]];
   Clear[bb0]; bb0 = {};
   dp = ((If[OptionValue[storeBB] == True, AppendTo[bb0, #]];
	  addTime[tdp, countdp++; dp0[#]])&);
   (* Compare sparse strategy against dense projection. *)
   If[OptionValue[debugProj],
     Module[{projDense = Nullspace[gauge]},
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
         Block[{y = dp[addTime[tdot, countdot++; hess.#]]},
               If[squareHess,
                  dp[addTime[tdot, countdot++; hess.y]],
                  y]]&,
         Length[grad],
	 Apply[Sequence, options],
	 (* Need to determine a reasonable estimate for this. *)
	 eigenPairs -> -Min[120, Length[grad]],
         orthoSubspace -> dp,
	 initialVector -> grad, printDetails -> 1];
     Module[{largeShiftSpace = cutoffNullspace[
         If[squareHess, Sqrt[vals], vals],
         vecs.grad, vecs, "zzz" -> zzz,
         Apply[Sequence, FilterRules[{opts}, Options[cutoffNullspace]]]]},
	    If[Length[largeShiftSpace] > 0,
	       (# - (largeShiftSpace.#).largeShiftSpace)&, Identity]]];
   smallProj = (addTime[tsp, countsp++; smallProj0[#]]&);
   (* Call the appropriate MINRES routine, setting some default values. *)
   result =
    minresFun[smallProj[dp[addTime[tdot, countdot++; hess.#]]]&,
     smallProj[dp[grad]], methodOptions[OptionValue[Method]],
     Apply[Sequence,
      FilterRules[{maxLanczosVecs -> 5000, printDetails -> 1},
       Options[minresFun]]]]; stepShifts = result[[1]];
   (* Verify result is in the subspace *)
   If[False,
      Print["Result overlap with small shift subspace, \
            dynamicPart, and both:  ",
	    {stepShifts.smallProj[stepShifts], stepShifts.dp[stepShifts],
	     stepShifts.smallProj[dp[stepShifts]]}/stepShifts.stepShifts]];
   Print["findDelta time (seconds):  dynamic part=", tdp,
         ", small shifts=", tsp, ", A=", tdot, ", total=",
         SessionTime[] - tinit];
   Print["findDelta counts:  dynamic part=", countdp,
         " (gauge transform=", dynamicPartCounter - countGauge,
         "), small shifts=", countsp, ", A=", countdot];
   If[OptionValue[storeHess],
      hess0 = hess; grad0 = grad; proj0 = Null; oo0 = Null;
      pairs0 = Null];
   -damping applyCutoff2[stepShifts, OptionValue[rescaleCutoff], zzz]] /;
 KeyExistsQ[minresLabels, methodName[OptionValue[Method]]];

SetAttributes[runRemote, HoldRest];
runRemote[command_, debug_, output_:None, error_:None] :=
    Block[{out},
          If[False, Print["Running:  ", command]];
          out = RunProcess[command,
               (* Path to scp, ssh ... *)
               ProcessEnvironment -> <|"PATH"->"/usr/bin:/usr/local/bin"|>];
          If[output =!= None,
             output = out["StandardOutput"]];
          If[debug > 0 && StringLength[out["StandardOutput"]]>0,
             Print[Style[out["StandardOutput"], FontColor -> Blue]]];
          If[error =!= None,
             error = out["StandardError"]];
          If[StringLength[out["StandardError"]]>0,
             Print[Style[out["StandardError"], FontColor -> Red]]];
          If[out["ExitCode"] != 0,
             Message[findDelta::external, out["ExitCode"]]];
          out["ExitCode"]];

(* External version.
   https://github.com/DaveGamble/cJSON
 *)
findDelta::poorMpiFlags = "These flag values are probably not very good.";
findDelta[data:{hess_, grad_, gauge_}, opts:OptionsPattern[]] :=
    Block[{
    damping = OptionValue[dampingFactor],
    action = OptionValue[externalAction],
    tinit = SessionTime[],
    symbolString = (a_Symbol -> b_) :> SymbolName[a] -> b,
    debug = printLevel[OptionValue[printDetails], 3],
    localPath = "",
    configFile = "hess-grad-gauge.json",
    shiftFile = "shifts.dat",
    outFile = "out.json",
    hessFile = "hess.mtx", gradFile = "grad.dat",
    gaugeFile = "gauge.mtx", out, remote,
    remoteHost = OptionValue[remoteHost],
    remotePath = "lattice/gauge-field-analyze/saddle-lib",
    executable = If[NumberQ[OptionValue[processes]],
                    "pshifts", "shifts"],
    mpi},
   remote = remoteHost <> ":" <> remotePath;
   mpi = If[NumberQ[OptionValue[processes]],
            Apply[Sequence,
                  {"mpirun",
                   "-np", ToString[OptionValue[processes]],
                   If[debug > 2, "--report-bindings", Nothing],
                   Which[OptionValue[mpiFlags] === None,
                         Nothing,
                         StringQ[OptionValue[mpiFlags]],
                         OptionValue[mpiFlags],
                         OptionValue[processes]<=4,
                         (* This seems to be the default for OpenMPI. *)
                         "--map-by numa", 
                         OptionValue[processes]<=8,
                         "--map-by l3cache",
                         (* In general, we want to bind np/8 consecutive
                           processes to each L3 cache.
                           The following does this for np=16. *)
                         True,
                         Message[findDelta::poorMpiFlags];
                         "--map-by core --bind-to l3cache"]}],
            Nothing];
   If[SymmetricMatrixQ[hess],
       Message[findDelta::symmetric]];
   If[action === "remote" || action == "detach" || action == "read",
      localPath = remoteHost <> "/";
      Run["mkdir -m 755 -p", localPath]];
   (* Dump dimensions, constants, and options into JSON file.

      Dump matrices and vectors: 
         hess, grad, gaugtransformationShifts
      The methods themselves will be chosen at compile time
      and we can always compare with the mathematica value.
    *)
   If[methodName[OptionValue[dynamicPartMethod]] =!= Automatic,
      Message[findDelta::dynamicPartMethod]; Return[$Failed]];
   If[action =!= "read",
      Export[localPath <> configFile, {
          "blockSize" -> nc^2 - 1,
          "partitions" -> If[OptionValue["partitions"] === Automatic,
                             latticeDimensions[[-1]],
                             OptionValue["partitions"]],
          "gaugePartitions" -> If[OptionValue["gaugePartitions"] === Automatic,
                                  Min[Length[gauge]/(nc^2 -1),
                                      latticeDimensions[[-1]]],
                                  OptionValue["gaugePartitions"]],
          "n" -> Length[grad],
          "gaugeDimension" -> Length[gauge],
          "removeConcave" -> OptionValue["removeConcave"],
          "eigenCutoffRescale" ->
              If[OptionValue[ignoreCutoff], 10.0^20, 1.0],
          "eigenCutoffMax" ->
              OptionValue["eigenCutoffMax"]/.Infinity -> 10.0^20,
          "eigenCutoff2" ->
              OptionValue["eigenCutoff2"]/.Infinity -> 10.0^20,
          "dynamicPartOptions" ->
              {methodOptions[OptionValue[dynamicPartMethod]]}/.symbolString,
          "largeShiftOptions" -> Append[
              OptionValue[largeShiftOptions]/.symbolString,
              "squareHess" -> !OptionValue["removeConcave"]],
          "linearSolveOptions" ->
	      {methodOptions[OptionValue[Method]]}/.symbolString,
          "chunkSize" -> OptionValue[chunkSize],
          "threads" -> OptionValue[threads],
          "printDetails" -> OptionValue[printDetails],
          "hessFile" -> hessFile,
          "gradFile" -> gradFile,
          "gaugeFile" -> gaugeFile}];
      (* If the matrix is symmetric, the default is to export
        only the lower triangle. *)
      Export[localPath <> hessFile, hess];
      Export[localPath <> gradFile, grad];
      Export[localPath <> gaugeFile, gauge];
      Apply[Clear, Unevaluated[data]]];
   (* Run external program locally *)
   If[action === Automatic,
      Run["rm", "-f", localPath <> shiftFile, localPath <> outFile];
      out = RunProcess[{mpi, "saddle-lib/"<>executable,
                        localPath <> configFile,
                        localPath <> shiftFile,
                        localPath <> outFile}];
      If[debug > 0,
         Print[Style[output = out["StandardOutput"], FontColor -> Blue]]];
      If[StringLength[out["StandardError"]]>0,
         Print[Style[error = out["StandardError"], FontColor -> Red]]];
      If[out["ExitCode"] != 0,
         Message[findDelta::external, out["ExitCode"]];
         Return[$Failed]]];
   (* Run external program on a remote host *)
   If[action == "remote" || action == "detach",
      If[debug > 2,
         Print["Sending files to ", remote]];
      If[runRemote[{"scp", "-B", localPath <> configFile,
                    localPath <> hessFile, localPath <> gradFile,
                    localPath <> gaugeFile,
                    remote}, debug] != 0, Return[$Failed]];
      If[runRemote[{"ssh", "-o", "BatchMode=true", remoteHost, "rm -f",
                    remotePath <> "/" <> shiftFile,
                    remotePath <> "/" <> outFile}, debug] != 0,
         Return[$Failed]]];
   (* Run external program on a remote host interactively *)
   If[action == "remote",
      If[runRemote[{"ssh", "-o", "BatchMode=true", remoteHost, mpi,
                    remotePath <> "/" <> executable,
                    remotePath <> "/" <> configFile,
                    remotePath <> "/" <> shiftFile,
                    remotePath <> "/" <> outFile},
                   debug, output, error] != 0, Return[$Failed]]];
   (* Start up external program and detach.
    After detaching, we then poll the external
    server every 60 seconds to see if the program is
    completed.  Thus, the connection to the remote server
    can be interrupted during the calculation. *)
   If[action == "detach",
      (* nohup is not relevant here, because there is no remote terminal.
        We just need to make sure stdout and stderr are redirected.
        For example:
            ssh localhost "(cd lattice/gauge-field-analyze/saddle-lib; touch junk; sleep 10s 1>/dev/null 2>/dev/null; rm junk) 1>/dev/null 2>&1 &"
       *)
      If[runRemote[{
          "ssh", "-o BatchMode=true", remoteHost,
          "(cd", remotePath, ";rm -f done;", mpi,
          " ./" <> executable, configFile, shiftFile, outFile,
          "1>shifts.log 2>shifts.err; touch done) >/dev/null 2>&1 &"
          }, debug] != 0, Return[$Failed]];
      If[debug>0, Print[
          DateObject[], "Started external program and detached."]]];
   If[action == "read" || action == "detach",
      Block[{dt = OptionValue[readInterval], tt = 0},
            While[RunProcess[{"ssh", "-o BatchMode=true", remoteHost,
                              "test -e",
                              remotePath <> "/done"},
                             "ExitCode"] != 0,
                  tt += dt;
                  Pause[dt]];
            If[debug>0, Print[
                DateObject[], "External program finished, ", tt, " s."]]];
      If[runRemote[{"scp",
                    remote <> "/" <> "shifts.log",
                    remote <> "/" <> "shifts.err",
                    localPath}, debug] != 0, Return[$Failed]];
      If[debug > 0,
         Print[Style[output = ReadString[localPath <> "shifts.log"],
                     FontColor -> Blue]]];
      error = ReadString[localPath <> "shifts.err"];
      If[error =!= EndOfFile,
         Print[Style[error, FontColor -> Red]]]];
   (* Read shifts from external file *)
   If[action == "read" || action == "remote" || action == "detach",
      If[runRemote[{"scp", remote <> "/" <> shiftFile,
                    remote <> "/" <> outFile,
                    localPath}, debug] != 0,
         Return[$Failed]]];
   stepOut = Import[localPath <> outFile];
   stepShifts = ReadList[localPath <> shiftFile, Number];
   -damping applyCutoff2[stepShifts, OptionValue[rescaleCutoff],
                         If[OptionValue[ignoreCutoff], Infinity, 1]]] /;
 methodName[OptionValue[Method]] === "External";

(*
  Use steepest descent to determine the next step.
  This assumes that grad is in the nullspace of gauge,
  so gauge is ignored.

  Eigenvalue-related cutoffs are ignored. 
 *)
findDelta::badGrad = "Gradient has overlap with gauge.";
findDelta[{hess_, grad_, gauge_}, OptionsPattern[]] :=
 Block[{zzz = If[OptionValue[ignoreCutoff], Infinity, 1],
        debug = printLevel[OptionValue[printDetails], 1],
        damping = OptionValue[dampingFactor], x},
    (* SymmetricMatrixQ[] has trouble with small nonzero elements *)
    If[matrixAsymmetry[hess] > 10^-12,
       Message[findDelta::asymmetric, matrixAsymmetry[hess]];
       Return[$Failed]];
    If[MatrixQ[gauge] &&
       Norm[gauge.grad] > 10^-12 Norm[grad] Norm[gauge],
       Message[findDelta::badGrad];
       Return[$Failed]];
    x = grad.grad/grad.hess.grad;
    stepShifts = x*grad;
    If[debug > 0,
       Print["steepest descent size:  ", x]];
    If[OptionValue[storeHess],
       Print["Define hess0, grad0"];
       hess0 = hess; grad0 = grad];
    -damping applyCutoff2[stepShifts, OptionValue[rescaleCutoff], zzz]] /;
 OptionValue[Method] === "steepestDescent";


SetAttributes[applyDelta, HoldFirst];
applyDelta[tgf_, delta_, fixed_] :=
  Block[{gradMap = gradLinkTake[fixed]},
   Do[
    Block[{coords = latticeCoordinates[i]},
     tgf[[dir, linearSiteIndex[coords]]] =
     getLink[tgf][dir, coords].
     If[boundaryLinkQ[dir, coords],
        IdentityMatrix[nc],
        MatrixExp[I Take[delta, gradMap[dir, coords]].SUGenerators[]]].
     getLink[tgf][dir, coords]],
    {dir, nd}, {i, latticeVolume[]}]];

latticeSaddlePointStep::usage =
  "Setting the option fixedDir=0 will use a projection onto the \
subspace orthogonal to all gauge transforms.  This should be the most \
correct method.  Can apply this to both dense and sparse systems.";
(* For passing options through to subfunctions: 
https://mathematica.stackexchange.com/questions/353/functions-with-options/ *)
Options[latticeSaddlePointStep] = Join[
    Options[findDelta], Options[actionHessian],
    {stepFile -> None, returnShifts -> False}];
latticeSaddlePointStep[opts:OptionsPattern[]] :=
 Block[{hess, grad, gauge, delta, gaugeField0,
        output = "", error = "",
        stepOut = None, stepShifts,
	t0 = SessionTime[], t1, t2, t3,
        debug = printLevel[OptionValue[printDetails], 3],
        action = OptionValue[externalAction]},
  gaugeField0 = makeRootLattice[];
  If[action =!= "read",
     {hess, grad} = actionHessian[
         getLink[gaugeField0],
         Apply[Sequence, FilterRules[{opts}, Options[actionHessian]]]];
     If[debug > 2,
        Print["Constructed hess, grad"]];
     gauge = If[OptionValue[fixedDir] > -1,
                gaugeTransformShifts[
                    gaugeField0, OptionValue[fixedDir]], None];
     If[debug > 2,
        Print["Constructed gauge"]]];
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
  delta = findDelta[
      {hess, grad, gauge},
      Apply[Sequence, FilterRules[{opts}, Options[findDelta]]]];
  (* Debug prints *)
  Which[
      False, Print["delta difference:  ", Chop[delta - delta0]],
      False, Print["delta:  ", delta],
      False, Print["delta, first link:  ", Take[delta, nc^2 - 1]]];
  t2 = SessionTime[];
  If[action =!= "write",
     applyDelta[gaugeField0, delta, OptionValue[fixedDir]];
     If[StringQ[OptionValue[stepFile]],
        If[debug > 2,
           Print["Saving step to ", OptionValue[stepFile]]];
        DeleteFile[OptionValue[stepFile]];
        Save[OptionValue[stepFile],
             {output, error, {opts}, stepOut,
              delta, stepShifts, gaugeField0}]]];
  t3 = SessionTime[];
  If[debug > 1,
     Print["latticeSaddlePointStep times:\n"
           <> "    matrices=", t1 - t0,
           " s, findDelta=", t2 - t1,
           " s, applyDelta=", t3 - t2, " s"]];
  If[OptionValue[returnShifts],
     {gaugeField0, stepShifts, stepOut},
     gaugeField0]];
