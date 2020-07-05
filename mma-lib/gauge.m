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

Options[landauGaugeHessian] = {fullMatrix -> False,
                               (* shift for calculating derivative *)
                               "delta" -> 10.0^-8};
landauGaugeHessian::usage = "Find the Hessian and gradient for the lattice \
norm (squared). \
Set the option fullMatrix to True to explicitly compute the upper \
triangle of the matrix.";
(*
  Method -> "simple" is a reference version with the code
  meant to be simple.
 *)
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
  Algebraic derivative
 *)
landauGaugeHessian[OptionsPattern[]] :=
 Block[{debug = True, count = 0},
  (* In the boundary case, make sure the lookup table is initialized
   before parallelization. *)
  If[latticeBC =!= "PERIODIC_GAUGEBC",
     reducedLinkIndex[1, Table[1, {nd}]]];
 (* Adding elements to a SparseArray one at a time is very inefficient
    in Mathematica; see
    https://mathematica.stackexchange.com/questions/777/efficient-by-element-updates-to-sparsearrays.
    Instead, we accumulate Array elements in an Association and create a
    SparseArray at the end. *)
 (*ParallelSum*) Sum[
   Block[{hess = Association[],
    grad = Array[0.0&, ngindex[]],
    gen = SUGenerators[] + 0.0 I,
    full = OptionValue[fullMatrix],
    idl = First[SUNorm[#]]^2&,
    coords},
   Do[coords = latticeCoordinates[i];
    Block[{
        uu = getLink[dir, coords],
        phases, vectors, center, vv, vvd, dd2, dd3,
        sl = coords2gindex[coords, 0],
        sr = coords2gindex[shift[dir, coords, -1], 0]
        },
     {phases, vectors, center} = getPhases[uu, True];
     vv = Transpose[vectors]; vvd = Conjugate[vectors];
     If[debug && False,
        Print["check vv ",
              Chop[uu - vv.DiagonalMatrix[Exp[I phases]].vvd]]];
     (* Use only the nc-1 diagonal generators *)
     dd2 = vv.Sum[gen[[j]]*(Diagonal[gen[[j]]].phases),
         {j, nc^2 -nc +1, nc^2 -1}].vvd;
     dd3 = 4*Map[reTrDot[#, dd2]&, gen];
     If[False,
        dd1 = Outer[If[Abs[#1 - #2] < 10^-16,
                       I, (#1-#2)/(Exp[I (#2 - #1)] - 1)]&, phases, phases];
        (* Should never access diagonal components *)
        Do[dd1[[i, i]] = Null, {i, nc}];
        Print["dd1 = ", MatrixForm[dd1]];
        dd2 = Table[If[i<nc^2 +1 -nc, gen[[i]]/dd1, I gen[[i]]], {i, nc^2-1}];
        Print["dd2 = ", ColumnForm[Map[MatrixForm, dd2]]]];
     If[debug && count++ < 1,
        Block[{delta = 10^-8, dl, dll, drr, dlr, zzz},
              dl = Table[(idl[MatrixExp[I delta gen[[ca1]]/2].uu] -
                             idl[MatrixExp[-I delta gen[[ca1]]/2].uu])/delta,
                            {ca1, nc^2 - 1}];
              zzz = 2 Table[Sum[Tr[gen[[i]].vvd.gen[[j]].vv]*
                                dl[[j]],
                                {j, nc^2 -1}],
                            {i, nc^2-1}];
              Print["Numeric gradient in diagonal space ", Chop[zzz]];
              Print["Numeric gradient ", dl];
              Print["algebraic gradient ", dd3];
              drr = Table[Sum[
                  z1*z2*idl[uu.MatrixExp[
                      -I delta (z1 gen[[ca1]] + z2 gen[[ca2]])/2]],
                  {z1, {-1, 1}}, {z2, {-1, 1}}]/delta^2,
                          {ca1, nc^2 - 1}, {ca2, nc^2 -1}];
              zzz = 4 Table[Sum[Tr[gen[[i1]].vvd.gen[[j1]].vv]*
                                Tr[gen[[i2]].vvd.gen[[j2]].vv]*
                                drr[[j1, j2]],
                                {j1, nc^2 -1}, {j2, nc^2 -1}],
                            {i1, nc^2-1}, {i2, nc^2-1}];
              Print["Numeric drr in diagonal space ",
                    MatrixForm[Chop[zzz]]];
              Print["Numeric drr ", MatrixForm[drr]];
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
              Print["Numeric dll in diagonal space ",
                    MatrixForm[Chop[zzz]]];
              Print["Numeric dll ", MatrixForm[dll]];
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
              Print["Numeric dlr in diagonal space ",
                    MatrixForm[Chop[zzz]]];
              Print["Numeric dlr ", MatrixForm[dlr]];
        ]];
     Do[
         (* Negative for inverse links. *)
         If[sr < Infinity,
            grad[[sl+ca1]] += dd3[[ca1]]];
         If[sl < Infinity,
            grad[[sr+ca1]] += -dd3[[ca1]]];
         t0 = SessionTime[];
         Do[
             oneAdd[hess, sr + ca1, sr + ca2,
                    Infinity,
                    full];
             oneAdd[hess, sl + ca1, sl + ca2,
                    0,
                    full];
             symAdd[hess, sl + ca1, sr + ca2,
                    0,
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
  {kernel, $KernelCount}]]/;OptionValue["delta"]==0;
