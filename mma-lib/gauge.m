(*
  Gauge transforms
 *)

Get[DirectoryName[$InputFileName] <> "matrix.m"];
coords2gindex[coords_, generator_] :=
  generator + (nc^2 - 1)*(
      If[boundarySiteQ[coords],
         Infinity,
         reducedSiteIndex[coords]] - 1);
ngindex[] := (nc^2 - 1)*reducedSiteCount[];

(* Gauge Transform Operator *)


(* Hessian and gradient for Landau gauge
  This is a modification of the function actionHessian
  in "saddle-point.m" *)

Options[landauGaugeHessian] = {
    fullMatrix -> False,
    (* delta=0:  calculate derivatives using algebraic expression.
      delta > 0: use finite differences to approximate derivatives.
      The two should match modulo numerical errors. *)
    "delta" -> 0,
    "debugHessian" -> False};
landauGaugeHessian::usage = "Find the Hessian and gradient for the lattice \
norm (squared). \
Set the option fullMatrix to True to explicitly compute the upper \
triangle of the matrix.";
landauGaugeHessian[OptionsPattern[]] :=
 Block[{delta = OptionValue["delta"]},
 (* Adding elements to a SparseArray one at a time is very inefficient
    in Mathematica; see
    https://mathematica.stackexchange.com/questions/777/efficient-by-element-updates-to-sparsearrays.
    Instead, we accumulate Array elements in an Association and create a
    SparseArray at the end. *)
  ParallelSum[
   Block[{hess = Association[],
    grad = Array[0.0&, ngindex[]],
    gen = SUGenerators[] + 0.0 I,
    full = OptionValue[fullMatrix],
    idl = First[SUNorm[#]]^2&,
    coords},
   Do[coords = latticeCoordinates[i];
    Block[{
        uu = getLink[dir, coords],
        sl = coords2gindex[coords, 0],
        sr = coords2gindex[shift[dir, coords, -1], 0]
        },
     Do[
         (* Negative for inverse links. *)
         If[sr < Infinity,
            grad[[sl+ca1]] +=
            (idl[MatrixExp[I delta gen[[ca1]]/2].uu] -
             idl[MatrixExp[-I delta gen[[ca1]]/2].uu])/delta];
         If[sl < Infinity,
            grad[[sr+ca1]] +=
            (idl[uu.MatrixExp[-I delta gen[[ca1]]/2]] -
             idl[uu.MatrixExp[I delta gen[[ca1]]/2]])/delta];
         t0 = SessionTime[];
         Do[
             oneAdd[hess, sr + ca1, sr + ca2,
                    Sum[
                        z1*z2*idl[uu.MatrixExp[
                            -I delta (z1 gen[[ca1]] + z2 gen[[ca2]])/2]],
                        {z1, {-1, 1}}, {z2, {-1, 1}}]/delta^2,
                    full];
             oneAdd[hess, sl + ca1, sl + ca2,
                    Sum[
                        z1*z2*idl[MatrixExp[
                            I delta (z1 gen[[ca1]] + z2 gen[[ca2]])/2].uu],
                        {z1, {-1, 1}}, {z2, {-1, 1}}]/delta^2,
                    full];
             symAdd[hess, sl + ca1, sr + ca2,
                    Sum[
                        z1*z2*idl[
                            MatrixExp[I delta z1 gen[[ca1]]/2].uu.
                                     MatrixExp[-I delta z2 gen[[ca2]]/2]],
                        {z1, {-1, 1}}, {z2, {-1, 1}}]/delta^2,
                    full],
	     {ca2, nc^2 - 1}],
         {ca1, nc^2 - 1}]],
      {dir, nd},
      {i, kernel, latticeVolume[], $KernelCount}];
   If[latticeBC =!= "PERIODIC_GAUGEBC",
      hess = KeySelect[hess,
                       (#[[1]]<Infinity && #[[2]]<Infinity)&]];
   {SparseArray[Normal[hess],
                {ngindex[], ngindex[]}], grad}],
  {kernel, $KernelCount}]]/;OptionValue["delta"]>0;

(*
  Algebraic derivatives
 *)
landauGaugeHessian[OptionsPattern[]] :=
    Block[{debug = OptionValue["debugHessian"], count = 0},
 (* Adding elements to a SparseArray one at a time is very inefficient
    in Mathematica; see
    https://mathematica.stackexchange.com/questions/777/efficient-by-element-updates-to-sparsearrays.
    Instead, we accumulate Array elements in an Association and create a
    SparseArray at the end. *)
  If[debug, Sum, ParallelSum][
   Block[{hess = Association[],
    grad = Array[0.0&, ngindex[]],
    gen = SUGenerators[] + 0.0 I,
    full = OptionValue[fullMatrix],
    coords},
   Do[coords = latticeCoordinates[i];
    Block[{
        uu = getLink[dir, coords],
        phases, vectors, center, vv, vvd, vvadj, grad2, grad3,
        adllrr, allrr, adlr, alr,
        sl = coords2gindex[coords, 0],
        sr = coords2gindex[shift[dir, coords, -1], 0]
        },
     {phases, vectors, center} = getPhases[uu, True];
     vv = Transpose[vectors]; vvd = Conjugate[vectors];
     (* Adjoint represenatation of vv.
      This code minimizes the number of matrix multiplications. *)
     vvadj = 2 Outer[reTrDot[#1, #2]&, Map[vv.#.vvd&, gen], gen, 1];
     If[debug && False,
        Print["vvadj = ", MatrixForm[vvadj]];
        Print["  verify orthogonality: ",
              MatrixForm[Chop[vvadj.Transpose[vvadj]]]]];
     If[debug && False,
        Print["check vv ",
              Chop[uu - vv.DiagonalMatrix[Exp[I phases]].vvd]]];
     (* This is just -I SULog[uu] *)
     grad2 = vv.(phases * vvd);
     grad3 = 2*Map[reTrDot[#, grad2]&, gen];
     (* left-left derivative in diagonalized link basis.
       Only construct the diagonal since off-diagonal elements are zero. *)
     adllrr = Table[
         (* Select non-diagonal vs. diagonal generators.
           In our ordering of the SU(N) generators,
           diagonal generators are at the end. *)
         If[i < nc^2 + 1 - nc,
            (* Only calculate one triangle, using symmetry. *)
            2*Sum[Block[
                {dphase = phases[[k1]]-phases[[k2]]},
                Re[gen[[i, k1, k2]] gen[[i, k2, k1]]]*dphase*
                   Cot[dphase/2]],
                  {k1, nc}, {k2, k1 - 1}],
            1],
         {i, nc^2 -1}];
     (* Equivalent to Transpose[vvadj].DiagonalMatrix[adllrr].vvadj *)
     allrr = LinearSolve[vvadj, adllrr*vvadj];
     (* left-right derivative in diagonalized link basis. *)
     adlr = Table[
         If[
             (* Select non-diagonal vs. diagonal generators. *)
             i < nc^2 + 1 - nc && j < nc^2 + 1 - nc,
             (* Our ordering of the SU(N) generators starts with all
               the real ones, then the imaginary ones.
               This picks out the non-zero cases. *)
             If[i ==j || Abs[i - j] == nc (nc-1)/2,
                (* Only calculate one triangle, using symmetry *)
                2*Sum[Block[
                    {dphase = phases[[k1]]-phases[[k2]]},
                    If[True,
                       (* Thee two are equivalent.
                        They differ by machine epsilon. *)
                       -dphase*If[i==j,
                                 Re[gen[[i, k1, k2]]*gen[[j, k2, k1]]]*
                                 Cot[dphase/2],
                                 Im[gen[[i, k1, k2]]*gen[[j, k2, k1]]]],
                       2*dphase*Im[gen[[i, k1, k2]]*gen[[j, k2, k1]]/
                                   (Exp[I dphase]-1)]]],
                      {k1, nc}, {k2, k1 - 1}],
                0],
             If[i == j, -1, 0]],
         {i, nc^2 -1}, {j, nc^2 - 1}];
     alr = LinearSolve[vvadj, adlr.vvadj];
     If[debug,
        Block[{delta = 10^-3, error = 10^-3, chopd = 10^-7,
               idl = First[SUNorm[#]]^2&,
               dderiv, ngrad, term1, term2,
               dl, dll, drr, dlr, zzz, yyy},
              dl = Table[(idl[MatrixExp[I delta gen[[ca1]]/2].uu] -
                             idl[MatrixExp[-I delta gen[[ca1]]/2].uu])/delta,
                            {ca1, nc^2 - 1}];
              zzz = 2 Table[Sum[Tr[gen[[i]].vvd.gen[[j]].vv]*
                                dl[[j]],
                                {j, nc^2 -1}],
                            {i, nc^2-1}];
              If[Max[Abs[dl-grad3]] > error,
                 Print["Numeric gradient in diagonal space ", Chop[zzz, chopd]];
                 Print["Numeric gradient ", dl];
                 Print["algebraic gradient ", grad3]];
              dderiv = Table[Sum[
                        z1*z2*SULog[MatrixExp[
                            I delta (z1 gen[[ca1]] + z2 gen[[ca2]])/2].
                                             DiagonalMatrix[Exp[I phases]]],
                        {z1, {-1, 1}}, {z2, {-1, 1}}]/delta^2,
                          {ca1, nc^2 - 1}, {ca2, nc^2 -1}];
              term2 = -2*I Map[Tr[#.DiagonalMatrix[phases]]&, dderiv, {2}];
              ngrad = Table[Sum[
                  z*SULog[MatrixExp[
                      I delta z gen[[ca]]/2].
                                   DiagonalMatrix[Exp[I phases]]],
                  {z, {-1, 1}}]/delta,
                             {ca, nc^2 - 1}];
              term1 = -2*Outer[Tr[#1.#2]&, ngrad, ngrad, 1];
              dll = Table[Sum[
                        z1*z2*idl[MatrixExp[
                            I delta (z1 gen[[ca1]] + z2 gen[[ca2]])/2].uu],
                        {z1, {-1, 1}}, {z2, {-1, 1}}]/delta^2,
                          {ca1, nc^2 - 1}, {ca2, nc^2 -1}];
              zzz = 4 Table[Sum[Tr[gen[[i1]].vvd.gen[[j1]].vv]*
                                Tr[gen[[i2]].vvd.gen[[j2]].vv]*
                                dll[[j1, j2]],
                                {j1, nc^2 -1}, {j2, nc^2 -1}],
                            {i1, nc^2-1}, {i2, nc^2-1}];
              If[Max[Abs[Flatten[allrr - dll]]] > error ||
                 (count++ < 1 && False),
                 Print["Second derivative term ",
                       MatrixForm[Chop[term2, chopd]]];
                 Print["First derivative term ",
                       MatrixForm[Chop[term1, chopd]]];
                 Print["Sum ",
                       MatrixForm[Chop[term1 + term2, chopd]]];
                 Print["Numeric dll in diagonal space ",
                       MatrixForm[Chop[zzz, chopd]]];
                 Print["algebraic dll in diagonal space",
                       MatrixForm[Chop[adllrr, chopd]]];
                 Print["Numeric dll ", MatrixForm[dll]];
                 Print["algebraic dll",
                       MatrixForm[Chop[allrr, chopd]]]];
              drr = Table[Sum[
                  z1*z2*idl[uu.MatrixExp[
                      -I delta (z1 gen[[ca1]] + z2 gen[[ca2]])/2]],
                  {z1, {-1, 1}}, {z2, {-1, 1}}]/delta^2,
                          {ca1, nc^2 - 1}, {ca2, nc^2 -1}];
              If[Max[Abs[Flatten[allrr - drr]]] > error,
                  Print["Numeric  drr ", MatrixForm[drr]]];
              dlr = Table[Sum[
                        z1*z2*idl[
                            MatrixExp[I delta z1 gen[[ca1]]/2].uu.
                                     MatrixExp[-I delta z2 gen[[ca2]]/2]],
                        {z1, {-1, 1}}, {z2, {-1, 1}}]/delta^2,
                          {ca1, nc^2 - 1}, {ca2, nc^2 -1}];
              zzz = 4 Table[Sum[Tr[gen[[i1]].vvd.gen[[j1]].vv]*
                                Tr[gen[[i2]].vvd.gen[[j2]].vv]*
                                dlr[[j1, j2]],
                                {j1, nc^2 -1}, {j2, nc^2 -1}],
                            {i1, nc^2-1}, {i2, nc^2-1}];
              If[Max[Abs[Flatten[alr - dlr]]] > error,
                 Print["Numeric dlr in diagonal space ",
                       MatrixForm[Chop[zzz, chopd]]];
                 Print["algebraic dlr in diagonal space",
                       MatrixForm[Chop[adlr, chopd]]];
                 Print["Numeric dlr ", MatrixForm[dlr]];
                 Print["algebraic dlr",
                       MatrixForm[Chop[alr, chopd]]]]
        ]];
     Do[
         (* Negative for inverse links. *)
         If[sr < Infinity,
            grad[[sl+ca1]] += grad3[[ca1]]];
         If[sl < Infinity,
            grad[[sr+ca1]] += -grad3[[ca1]]];
         t0 = SessionTime[];
         Do[
             oneAdd[hess, sr + ca1, sr + ca2,
                    allrr[[ca1, ca2]],
                    full];
             oneAdd[hess, sl + ca1, sl + ca2,
                    allrr[[ca1, ca2]],
                    full];
             symAdd[hess, sl + ca1, sr + ca2,
                    alr[[ca1, ca2]],
                    full],
	     {ca2, nc^2 - 1}],
         {ca1, nc^2 - 1}]],
      {dir, nd},
      {i, kernel, latticeVolume[], If[debug, 1, $KernelCount]}];
   If[latticeBC =!= "PERIODIC_GAUGEBC",
      hess = KeySelect[hess,
                       (#[[1]]<Infinity && #[[2]]<Infinity)&]];
   {SparseArray[Normal[hess],
                {ngindex[], ngindex[]}], grad}],
   {kernel, If[debug, 1, $KernelCount]}]]/;OptionValue["delta"]==0;

(*
  Create matrix associated with global color rotations
 *)
globalColorRotations::usage =
  "Return a matrix generating the nc^2-1 global color rotations.  \
The result is returned as a SparseArray.";
globalColorRotations[] :=
    SparseArray[
        Table[
            {Mod[i-1, nc^2 -1]+1, i} -> 1,
            {i, ngindex[]}],
        {nc^2-1, ngindex[]}];

SetAttributes[applyGaugeTransform, HoldFirst];
applyGaugeTransform[gf_, delta_] :=
 Block[{gi, gg},
  Do[
   Block[{coords = latticeCoordinates[i]},
    gi = coords2gindex[coords, 0] + {1, nc^2 -1};
    gg = MatrixExp[I Take[delta, gi].SUGenerators[]];
    Do[
        gf[[dir, linearSiteIndex[coords]]] =
        gg.getLink[gf][dir, coords];
        gf[[dir, linearSiteIndex[shift[dir, coords, -1]]]] =
        getLink[gf][dir, shift[dir, coords, -1]].ConjugateTranspose[gg],
        {dir, nd}]],
   {i, latticeVolume[]}]];

minimumNormGaugeStep::usage =
  "Find a gauge transform that minimizes the lattice norm.";
(* This is a modification of latticeSaddlePointStep. *)
Options[minimumNormGaugeStep] = Join[
    Options[findDelta], Options[landauGaugeHessian],
    {stepFile -> None, returnShifts -> False}];
minimumNormGaugeStep[opts:OptionsPattern[]] :=
 Block[{hess, grad, gauge, delta,
        output = "", error = "",
        stepOut = None, stepShifts,
	t0 = SessionTime[], t1, t2, t3,
        debug = printLevel[OptionValue[printDetails], 3],
        action = OptionValue[externalAction]},
  gaugeField0 = gaugeField;
  If[action =!= "read",
     Print["calling landauGaugeHessian"];
     {hess, grad} = landauGaugeHessian[
         Apply[Sequence, FilterRules[{opts}, Options[landauGaugeHessian]]]];
     If[debug > 2,
        Print["Constructed hess, grad"]];
     gauge = globalColorRotations[];
     If[debug > 2,
        Print["Constructed gauge"]]];
  t1 = SessionTime[];
  (* Debug prints *)
  Which[False,
        Print["hessian difference:  ", Chop[Normal[hess] - hess0]];
        Print["gradient difference:  ", Chop[grad - gradient0]],
        False,
        Print["hessian:", Normal[hess]];
        Print["gradient:  ", grad];
        Print["gauge:  ", gauge],
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
     applyGaugeTransform[gaugeField0, delta];
     If[StringQ[OptionValue[stepFile]],
        If[debug > 2,
           Print["Saving step to ", OptionValue[stepFile]]];
        DeleteFile[OptionValue[stepFile]];
        Save[OptionValue[stepFile],
             {output, error, {opts}, stepOut,
              delta, stepShifts, gaugeField0}]]];
  t3 = SessionTime[];
  If[debug > 1,
     Print["minimumNormGaugeStep times:\n"
           <> "    matrices = ", t1 - t0,
           " s, findDelta = ", t2 - t1,
           " s, applyGaugeTransform = ", t3 - t2, " s"]];
  If[OptionValue[returnShifts],
     {gaugeField0, stepShifts, stepOut},
     gaugeField0]];
