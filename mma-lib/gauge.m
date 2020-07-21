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
    (* Set to false to compare with single-site minimization *)
    "includeLeftRight" -> True,
    "timing" -> False, "debugHessian" -> False};
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
             If[OptionValue["includeLeftRight"],
                symAdd[hess, sl + ca1, sr + ca2,
                    Sum[
                        z1*z2*idl[
                            MatrixExp[I delta z1 gen[[ca1]]/2].uu.
                                     MatrixExp[-I delta z2 gen[[ca2]]/2]],
                        {z1, {-1, 1}}, {z2, {-1, 1}}]/delta^2,
                       full]],
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
  Algebraic derivatives, debug version
 *)
landauGaugeHessian[OptionsPattern[]] :=
 Block[{
     hess = Association[],
     grad = Array[0.0&, ngindex[]],
     gen = SUGenerators[],
     full = OptionValue[fullMatrix],
     diagLookup, offDiagLookup,
     coords, count = 0},
  (* Non-diagonal generators of SU(N) *)
  diagLookup = Flatten[Table[
      If[gen[[i, k1, k2]] != 0 && k1 != k2,
         (* factor of 2 since we only do k2 < k2 *)
         {i, k1, k2, N[2*gen[[i, k1, k2]] gen[[i, k2, k1]]]},
         Nothing],
      {i, nc^2-1}, {k1, 2, nc}, {k2, k1 -1}], 2];
  offDiagLookup = Flatten[Table[
      If[Im[gen[[i, k1, k2]] gen[[j, k2, k1]]] != 0 &&
         k1 != k2,
         {i, j, k1, k2, N[2*Im[gen[[i, k1, k2]] gen[[j, k2, k1]]]]},
         Nothing],
      {i, 2, nc^2-1}, {j, i-1}, {k1, 2, nc}, {k2, k1 -1}], 3];
  Do[coords = latticeCoordinates[i];
   Block[{
     uu = getLink[dir, coords],
     phases, vectors, center, vv, vvd, vvadj, grad2, grad3,
     adllrr, allrr, adlr, alr,
     sl = coords2gindex[coords, 0],
     sr = coords2gindex[shift[dir, coords, -1], 0]},
    {phases, vectors, center} = getPhases[uu, True];
    vv = Transpose[vectors]; vvd = Conjugate[vectors];
    (* Adjoint represenatation of vv.
      This code minimizes the number of matrix multiplications. *)
    vvadj = 2 Outer[reTrDot[#1, #2]&, Map[vv.#.vvd&, gen], gen, 1];
    If[False,
       Print["vvadj = ", MatrixForm[vvadj]];
       Print["  verify orthogonality: ",
             MatrixForm[Chop[vvadj.Transpose[vvadj]]]]];
    If[False,
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
        (count++ < 1 && True),
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
     If[Max[Abs[Flatten[alr - dlr]]] > error || count<2,
        Print["Numeric dlr in diagonal space ",
              MatrixForm[Chop[zzz, chopd]]];
        Print["algebraic dlr in diagonal space",
              MatrixForm[Chop[adlr, chopd]]];
        Print["Numeric dlr ", MatrixForm[dlr]];
        Print["algebraic dlr",
              MatrixForm[Chop[alr, chopd]]]]
    ];
    Do[
        (* Negative for inverse links. *)
        If[sr < Infinity,
           grad[[sl+ca1]] += grad3[[ca1]]];
        If[sl < Infinity,
           grad[[sr+ca1]] += -grad3[[ca1]]];
        Do[
            oneAdd[hess, sr + ca1, sr + ca2,
                   allrr[[ca1, ca2]],
                   full];
            oneAdd[hess, sl + ca1, sl + ca2,
                   allrr[[ca1, ca2]],
                   full];
            If[OptionValue["includeLeftRight"],
               symAdd[hess, sl + ca1, sr + ca2,
                   alr[[ca1, ca2]],
                      full]],
	    {ca2, nc^2 - 1}],
        {ca1, nc^2 - 1}]],
     {dir, nd},
     {i, latticeVolume[]}];
  If[latticeBC =!= "PERIODIC_GAUGEBC",
     hess = KeySelect[hess,
                      (#[[1]]<Infinity && #[[2]]<Infinity)&]];
  {SparseArray[Normal[hess],
               {ngindex[], ngindex[]}], grad}
 ]/; OptionValue["delta"]==0 && OptionValue["debugHessian"];


(*
  Algebraic derivatives, fast
 *)
(* Adding elements to a SparseArray one at a time is very inefficient
  in Mathematica; see
  https://mathematica.stackexchange.com/questions/777/efficient-by-element-updates-to-sparsearrays.
  Instead, we accumulate Array elements in an Association and create a
  SparseArray at the end. *)
landauGaugeHessian[OptionsPattern[]] :=
   ParallelSum[
   Block[{hess = Association[],
     grad = Array[0.0&, ngindex[]],
     gen = SUGenerators[],
     full = OptionValue[fullMatrix],
     includeLeftRight = OptionValue["includeLeftRight"],
     diagLookup, offDiagLookup,
          addTimeNull = If[
              OptionValue["timing"],
              Function[{timer, expr},
                       Block[{t1 = SessionTime[]},
                             expr; timer += SessionTime[] - t1],
                       {HoldAll, SequenceHold}],
              Null&],
     t0 = SessionTime[], tm = 0, ta = 0, tc = 0,
     coords},
    diagLookup = Table[Block[
        (* For Cartan (diagonal) generators. *)
        {z = {0, 0, 2 Tr[gen[[i]].gen[[i]]]}},
        Do[
            If[gen[[i, k1, k2]] != 0,
               (* factor of 2 since we only do k1 < k2 *)
               z = {k1, k2, N[Re[2*gen[[i, k1, k2]]*gen[[i, k2, k1]]]]}],
            {k1, 2, nc}, {k2, k1 -1}];
        z], {i, nc^2-1}];
   offDiagLookup = Flatten[Table[
       If[Im[gen[[i, k1, k2]] gen[[j, k2, k1]]] != 0,
          (* factor of 2 since we only do k1 < k2 *)
          {i, j, k1, k2, N[-2*Im[gen[[i, k1, k2]] gen[[j, k2, k1]]]]},
          Nothing],
       {i, nc^2-1}, {j, nc^2 -1}, {k1, 2, nc}, {k2, k1 -1}], 3];
   Do[coords = latticeCoordinates[i];
    Block[{
        uu = getLink[dir, coords],
        phases, vectors, center,
        vv, vvd, vvadj, grad2, grad3,
        adllrr, allrr, adlr, alr,
        sl = coords2gindex[coords, 0],
        sr = coords2gindex[shift[dir, coords, -1], 0]
        },
     {phases, vectors, center} = getPhases[uu, True];
     vv = Transpose[vectors]; vvd = Conjugate[vectors];
     (* Adjoint represenatation of vv.
      This code minimizes the number of matrix multiplications. *)
     vvadj = 2 Outer[reTrDot[#1, #2]&, Map[vv.#.vvd&, gen], gen, 1];
     (* This is just -I SULog[uu] *)
     grad2 = vv.(phases * vvd);
     grad3 = 2*Map[reTrDot[#, grad2]&, gen];
     (* left-left derivative in diagonalized link basis.
       Only construct the diagonal since off-diagonal elements are zero. *)
     addTimeNull[tc, adllrr = Map[
         If[#[[1]] == 0,
            #[[3]],
            Block[{dphase = phases[[#[[1]]]]-phases[[#[[2]]]]},
                  #[[3]]*dphase*Cot[dphase/2]]]&,
           diagLookup]];
     (* Equivalent to Transpose[vvadj].DiagonalMatrix[adllrr].vvadj *)
     addTimeNull[tm, allrr = LinearSolve[vvadj, adllrr*vvadj]];
     (* left-right derivative in diagonalized link basis. *)
     attTimeNull[tc, adlr = SparseArray[Map[
         (Take[#, 2] -> #[[5]]*(phases[[#[[3]]]]-phases[[#[[4]]]]))&,
         offDiagLookup], {nc^2 - 1, nc^2 -1}]];
     (* Equivalent to
       Transpose[vvadj].(adlr - DiagonalMatrix[adllrr]).vvadj *)
     addTimeNull[tm, alr = LinearSolve[vvadj, adlr.vvadj - adllrr*vvadj]];
     addTimeNull[ta, Do[
         (* Negative for inverse links. *)
         If[sr < Infinity,
            grad[[sl+ca1]] += grad3[[ca1]]];
         If[sl < Infinity,
            grad[[sr+ca1]] += -grad3[[ca1]]];
         Do[
             oneAdd[hess, sr + ca1, sr + ca2,
                    allrr[[ca1, ca2]],
                    full];
             oneAdd[hess, sl + ca1, sl + ca2,
                    allrr[[ca1, ca2]],
                    full];
             If[includeLeftRight, symAdd[hess, sl + ca1, sr + ca2,
                    alr[[ca1, ca2]],
                    full]],
	     {ca2, nc^2 - 1}],
         {ca1, nc^2 - 1}]]],
      {dir, nd},
      {i, kernel, latticeVolume[], $KernelCount}];
   If[OptionValue["timing"],
      Print[kernel, ": times: ", {tc, tm, ta, SessionTime[] - t0}]];
   If[latticeBC =!= "PERIODIC_GAUGEBC",
      hess = KeySelect[hess,
                       (#[[1]]<Infinity && #[[2]]<Infinity)&]];
   {SparseArray[Normal[hess],
                {ngindex[], ngindex[]}], grad}],
   {kernel, $KernelCount}]/;
    OptionValue["delta"]==0 && !OptionValue["debugHessian"];


(*
  Create matrix associated with global color rotations
 *)
globalColorRotations::usage =
  "Return a matrix generating the nc^2-1 global color rotations.  \
The result is returned as a SparseArray.";
globalColorRotations[] := Block[{z = N[1/Sqrt[latticeVolume[]]]},
    SparseArray[
        Table[
             (* saddle-lib expects floating point matrix elements *)
            {Mod[i-1, nc^2 -1]+1, i} -> z,
            {i, ngindex[]}],
        {nc^2-1, ngindex[]}]];

SetAttributes[applyGaugeTransform, HoldFirst];
applyGaugeTransform[gf_, delta_] :=
 Block[{gi, gg, sum, asum, count=0},
  Do[
   Block[{coords = latticeCoordinates[i]},
    gi = coords2gindex[coords, 0] + {1, nc^2 -1};
    gg = MatrixExp[I Take[delta, gi].SUGenerators[]];
    If[False, 
       sum = Sum[
           SULog[getLink[dir, coords]]
           - SULog[getLink[dir, shift[dir, coords, -1]]],
           {dir, nd}];
       asum = -Map[2*Im[Tr[sum.#]]&, SUGenerators[]]/(2 nd);
       If[count++ < 10,
          Print["norms: ",
             {Norm[Take[delta, gi] - asum],
              Take[delta, gi].asum/(Norm[Take[delta, gi]]*Norm[asum]),
              Norm[Take[delta, gi]],
              Norm[asum]}]]];
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
        action = OptionValue[externalAction],
        gaugeField0 = gaugeField},
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
