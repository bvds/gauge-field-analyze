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
    "debug" -> False};
landauGaugeHessian::usage = "Find the Hessian and gradient for the lattice \
norm (squared). \
Set the option fullMatrix to True to explicitly compute the upper \
triangle of the matrix.";
landauGaugeHessian[OptionsPattern[]] :=
 Block[{delta = OptionValue["delta"]},
  (* In the boundary case, make sure the lookup table is initialized
   before parallelization. *)
  If[latticeBC =!= "PERIODIC_GAUGEBC",
     reducedLinkIndex[1, Table[1, {nd}]]];
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
    Block[{debug = OptionValue["debug"], count = 0},
  (* In the boundary case, make sure the lookup table is initialized
   before parallelization. *)
  If[latticeBC =!= "PERIODIC_GAUGEBC",
     reducedLinkIndex[1, Table[1, {nd}]]];
 (* Adding elements to a SparseArray one at a time is very inefficient
    in Mathematica; see
    https://mathematica.stackexchange.com/questions/777/efficient-by-element-updates-to-sparsearrays.
    Instead, we accumulate Array elements in an Association and create a
    SparseArray at the end. *)
  If[debug, Sum, ParallelSum][
   Block[{hess = Association[],
    grad = Array[0.0&, ngindex[]],
    gen = SUGenerators[] + 0.0 I,
    full = OptionValue[fullMatrix]
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
     (* left-left derivative in diagonalized link basis. *)
     adllrr = Table[
         (* Select non-diagonal vs. diagonal generators.
           In our ordering of the SU(N) generators,
           diagonal generators are at the end. *)
         If[i < nc^2 + 1 - nc,
            Sum[Block[{dphase = phases[[k1]]-phases[[k2]]},
                      Abs[gen[[i, k1, k2]]] dphase*
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
             (* Specific to our ordering of the SU(N) generators:
               First all the real ones, then the imaginary. *)
             If[i ==j || Abs[i - j] == nc (nc-1)/2,
                Sum[Block[{dphase = phases[[k1]]-phases[[k2]]},
                          Abs[gen[[i, k1, k2]]]*dphase*
                             If[i==j, -Cot[dphase/2], Sign[j - i]]],
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
