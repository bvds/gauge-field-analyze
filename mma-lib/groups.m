(* Many of these functions use the global variable "nc" *)
(* One can't use the pattern n_:nc because nc may change. *)

randomSUMatrix[] := randomSUMatrix[nc];
randomSUMatrix[n_] :=
 Block[{mat = RandomVariate[CircularUnitaryMatrixDistribution[n]]},
       mat*Exp[2 Pi I RandomInteger[n-1]/n]/Det[mat]^(1/n)];
centerSUMatrix[j_] := centerSUMatrix[j, nc];
centerSUMatrix[j_, n_] := IdentityMatrix[n] Exp[2 Pi I j/n];
zeroMatrix[] := zeroMatrix[nc];
zeroMatrix[n_] := Array[0&, {n, n}];
Options[matrixEqual] = {Tolerance -> $MachineEpsilon, "prefix" -> None};
matrixEqual[m1_, m2_, OptionsPattern[]] :=
 Block[{norm = Norm[Flatten[m1 - m2], Infinity]},
       If[norm < OptionValue[Tolerance] Sqrt[Length[m1] Length[m1[[1]]]],
          True,
          Print[If[StringQ[OptionValue["prefix"]],
                   OptionValue["prefix"] <> ":  ", ""] <>
                "Matrices not equal, norm:  ", norm]; False]];
Options[SUMatrixQ] = Options[matrixEqual];
SUMatrixQ[mat_?SquareMatrixQ, OptionsPattern[]] :=
 Block[{diffs = {Det[mat] - 1,
     mat.ConjugateTranspose[mat] - IdentityMatrix[Length[mat]]},
   norm}, norm = Norm[Flatten[diffs], Infinity];
          If[norm < OptionValue[Tolerance] Length[mat], True,
             Print[If[StringQ[OptionValue["prefix"]],
                      OptionValue["prefix"] <> ":  ", ""] <>
                     "Not SU matrix, error:  ", norm];
             Print["Details:  ", diffs];
             False]];

SUGenerators::usage = "Construct the generators for SU(N), caching the result.
Normalization Tr[T_a T_b] = delta_{ab}/2.";
SUGenerators[] := SUGenerators[nc];
SUGenerators[n_?IntegerQ] :=
 SUGenerators[n] =
  Join[Flatten[
    Table[If[(i == k1 && j == k2) || (i == k2 && j == k1), 1/2,
      0], {i, 2, n}, {j, i - 1}, {k1, n}, {k2, n}], 1],
   Flatten[Table[
     Which[i == k1 && j == k2, I/2, i == k2 && j == k1, -I/2, True,
      0], {i, 2, n}, {j, i - 1}, {k1, n}, {k2, n}], 1],
   Table[Which[k1 != k2, 0, k1 < i, 1/Sqrt[2 i (i - 1)],
     k1 == i, (1 - i)/Sqrt[2 i (i - 1)], True, 0], {i, 2, n}, {k1,
	n}, {k2, n}]];
SUSymmetric::usage = "Construct the symmetric coefficients d_{a,b,c},
caching the result.  This tensor is sparse, but we will be inefficient for now.";
SUSymmetric[] := SUSymmetric[nc];
SUSymmetric[n_?IntegerQ] :=
 SUSymmetric[n] =
  Block[{gen = SUGenerators[n]},
	2 Outer[Tr[(#1.#2 + #2.#1).#3]&, gen, gen, gen, 1]];

centerPhases[mat_] :=(* Used for some testing *)
 Block[{delta, phases = Arg[Eigenvalues[mat]]},
  phases[[1]] -= Apply[Plus, phases];
  Do[If[phases[[i]] - phases[[j]] > 2*Pi,
    delta = Floor[(phases[[i]] - phases[[j]] + 2*Pi)/(4*Pi)];
    phases[[i]] -= 2*Pi*delta; phases[[j]] += 2*Pi*delta], {i,
    Length[phases]}, {j, Length[phases]}];
  Sort[phases, Order[Abs[#1], Abs[#2]]&]];

Options[getPhases] = {Tolerance -> 10^-5, "debug"->False};
getPhases::phaseSum = "Invalid determinant, total phase = `1`";
getPhases[mat_, OptionsPattern[]]:= Block[{values, vectors, phases, phaseSum,
                         fix=True, delta},
  {values, vectors} = Eigensystem[mat];
  phases = Arg[values]; phaseSum = Total[phases];
  (* Gram-Schmidt orthogonalization, for algebraic cases. *)
  vectors = Orthogonalize[vectors];
  (* Sanity check.  QDP++ default is single precision,
  so we set the default tolerance high. *)
  If[Abs[Mod[phaseSum + Pi, 2 Pi] - Pi] > OptionValue[Tolerance],
   Message[getPhases::phaseSum, phaseSum]; Return[$Failed]];
  (* If the total is nonzero, need to shift *)
  If[phaseSum != 0, phases[[1]] -= phaseSum;
   While[fix, fix = False;
    Do[If[phases[[i]] - phases[[j]] > 2*Pi,
      delta = Floor[(phases[[i]] - phases[[j]] + 2*Pi)/(4*Pi)];
      phases[[i]] -= 2*Pi*delta; phases[[j]] += 2*Pi*delta;
      If[OptionValue["debug"], Print["Fixing phases ", {i,j}]];
      fix = True], {i, Length[mat]}, {j, Length[mat]}]]];
  {phases, vectors}];
SUPower::usage = "Take some power of an SU(N) matrix with the property that U^z, as a function of z is a smooth map onto the group manifold.  The analogous Mathematica function MatrixPower[] uses a different convention for the branch taken in the complex plane.";
Options[SUPower] = Options[getPhases];
SUPower[mat_, power_, opts:OptionsPattern[]] :=
 Block[{vectors, phases, result},
  {phases, vectors} = getPhases[mat, opts];
  result = Transpose[vectors].DiagonalMatrix[
     Exp[I*phases*power]].Conjugate[vectors];
  If[OptionValue["debug"] && Not[SUMatrixQ[result, OptionValue[Tolerance]]],
     Print["matrix power not SU"]];
  result];
SULog::usage = "Take the logaritm of an SU(N) matrix, choosing the root such that the resulting matrix is traceless and has eigenvalues that obey Abs[lamda_i-lambda_j] <= 2 Pi.  The analogous Mathematica function MatrixLog[] uses a different convention for the branch taken in the complex plane.";
Options[SULog] = Options[getPhases];
SULog[mat_, opts:OptionsPattern[]] :=
 Block[{vectors, phases},
  {phases, vectors} = getPhases[mat, opts];
  Transpose[vectors].DiagonalMatrix[I phases].Conjugate[vectors]];
