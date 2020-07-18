(* Many of these functions use the global variable "nc" *)
(* One can't use the pattern n_:nc because nc may change. *)

randomSUMatrix[] := randomSUMatrix[nc];
randomSUMatrix[n_] :=
 Block[{mat = RandomVariate[CircularUnitaryMatrixDistribution[n]]},
  (* Arg[Det[mat]] is uniformly distributed on [-Pi, Pi],
    so the distribution over SU(N) is uniform. *)
       mat/Det[mat]^(1/n)];
randomPhaseDistribution[x_List] := 
  Block[{nc = Length[x]}, 
   2^((nc - 1)*(nc - 2)/2)/Pi^(nc - 1)*
   Product[1 - Cos[x[[i]] - x[[j]]], {i, nc - 1}, {j, i + 1, nc}]];
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

trace::usage = "Trace of the product of Hermitian matrices.  \
This allows for future optimization.";
trace[a_, b_] := Block[{x = Flatten[a].Flatten[Transpose[b]]},
                       If[Abs[Im[x]]>10^-8,
                          Print["non-hermetian arg", {a, b}];
                          Print[Column[Map[Short, Stack[_]]]];
                          Abort[]];
                       Re[x]];
trace[a__] := Tr[Dot[a]]/;Length[{a}]!=2;

SUGenerators::usage = "Construct the generators for SU(N), caching the result.
Normalization Tr[T_a T_b] = delta_{ab}/2.";
SUGenerators[] := SUGenerators[nc];
SUGenerators[n_?IntegerQ] := SUGenerators[n] =
Join[
    Flatten[Table[
        If[(i == k1 && j == k2) || (i == k2 && j == k1), 1/2, 0],
        {i, 2, n}, {j, i - 1}, {k1, n}, {k2, n}], 1],
    Flatten[Table[
        Which[i == k1 && j == k2, I/2, i == k2 && j == k1, -I/2, True, 0],
        {i, 2, n}, {j, i - 1}, {k1, n}, {k2, n}], 1],
    Table[Which[k1 != k2, 0, k1 < i, 1/Sqrt[2 i (i - 1)],
                k1 == i, (1 - i)/Sqrt[2 i (i - 1)], True, 0],
          {i, 2, n}, {k1, n}, {k2, n}]];
SUSymmetric::usage = "Construct the symmetric coefficients d_{a,b,c},
caching the result.  This tensor is sparse.";
SUSymmetric[] := SUSymmetric[nc];
SUSymmetric[n_?IntegerQ] :=
 SUSymmetric[n] =
  Block[{gen = SUGenerators[n]},
	2 SparseArray[Outer[Tr[(#1.#2 + #2.#1).#3]&, gen, gen, gen, 1]]];
SUStructure::usage = "Construct the structure constants f_{a,b,c},
caching the result.  This tensor is sparse.";
SUStructure[] := SUStructure[nc];
SUStructure[n_?IntegerQ] :=
 SUStructure[n] =
  Block[{gen = SUGenerators[n]},
	2 SparseArray[Im[Outer[Tr[(#1.#2 - #2.#1).#3]&, gen, gen, gen, 1]]]];

centerPhases[mat_] :=(* Used for some testing *)
 Block[{delta, phases = Arg[Eigenvalues[mat]]},
  phases[[1]] -= Apply[Plus, phases];
  Do[If[phases[[i]] - phases[[j]] > 2*Pi,
    delta = Floor[(phases[[i]] - phases[[j]] + 2*Pi)/(4*Pi)];
    phases[[i]] -= 2*Pi*delta; phases[[j]] += 2*Pi*delta],
     {i, Length[phases]}, {j, Length[phases]}];
  Sort[phases, Greater]];

equalPhase::usage = "Determine whether two sets of phases represent the same element of the group.  Mostly for debugging:  verify that cleanPhases[] and getPhases[] handle the boundary cases correctly.";
equalPhase[x_,y_] := MatrixExp[DiagonalMatrix[I x]] ==
                     MatrixExp[DiagonalMatrix[I y]];

cleanPhases[phases0_] :=
 Block[{phases = phases0, n = Length[phases0], fix = True},
       phases[[1]] -= Total[phases];
       While[
           fix,
           fix = False;
           Do[If[
               phases[[i]] - phases[[j]] > 2*Pi,
               delta = Floor[(phases[[i]] - phases[[j]] + 2*Pi)/(4*Pi)];
               phases[[i]] -= 2*Pi*delta; phases[[j]] += 2*Pi*delta;
               fix = True],
              {i, n}, {j, n}]];
       phases];

getPhases::phaseSum = "Invalid determinant, total phase = `1`";
Options[getPhases] = {Tolerance -> 10^-5, "center" -> False, "debug" -> False};
(* Create a class-like structure with private variables.
  This is a work-around for the fact that compiled functions
  cannot return multiple variable types. *)
Module[
  {fixPhases, bestPhases, bestCenter},
 fixPhases = Compile[
  {{values, _Complex, 1}, {tolerance, _Real}, {centerFlag, True|False}},
   Block[{nc = Length[values], phases, phaseSum, fix, delta},
    (* If the total is nonzero, need to shift *)
    phases = Arg[values];
    phaseSum = Total[phases];
    (* Sanity check.  QDP++ default is single precision,
      so we set the default tolerance high. *)
    If[Abs[Mod[phaseSum + Pi, 2 Pi] - Pi] > tolerance,
       Message[getPhases::phaseSum, phaseSum]];
    Do[
        If[center > 0, phases -= 2*Pi/nc; phaseSum = -2 Pi];
        phases[[1]] -= phaseSum;
        fix = True;
        While[
            fix,
            fix = False;
            Do[If[
                phases[[i]] - phases[[j]] > 2*Pi,
                delta = Floor[(phases[[i]] - phases[[j]] + 2*Pi)/(4*Pi)];
                phases[[i]] -= 2*Pi*delta; phases[[j]] += 2*Pi*delta;
                fix = True],
               {i, nc}, {j, nc}]];
        If[center == 0 || Norm[bestPhases] > Norm[phases],
           bestPhases = phases; bestCenter = center],
        {center, 0, If[centerFlag, nc-1, 0]}]]];
 getPhases[mat_, matrixFlag_, OptionsPattern[]] :=
    Block[{values, vectors},
          If[matrixFlag,    
             {values, vectors} = Eigensystem[mat];
             (* Gram-Schmidt orthogonalization, for algebraic cases. *)
             vectors = Orthogonalize[vectors],
             values = Eigenvalues[mat];
             vectors = Nothing];
          fixPhases[values, OptionValue[Tolerance], OptionValue["center"]];
          {bestPhases, vectors, bestCenter}]];
SUPower::usage = "Take some power of an SU(N) matrix with the property that U^z, as a function of z is a smooth map onto the group manifold.  The analogous Mathematica function MatrixPower[] uses a different convention for the branch taken in the complex plane.";
Options[SUPower] = Options[getPhases];
SUPower[mat_?MatrixQ, power_, opts:OptionsPattern[]] :=
 Block[{vectors, phases, result, center},
  {phases, vectors, center} = getPhases[mat, True, opts];
  result = Transpose[vectors].(
      (* equivalent to Diagonal[Exp[I phases]].Conjugate[vectors] *)
      Exp[I*phases*power]*Conjugate[vectors]);
  If[OptionValue["debug"] && Not[SUMatrixQ[result, OptionValue[Tolerance]]],
     Print["matrix power not SU"]];
  result];
SULog::usage = "Take the logaritm of an SU(N) matrix, choosing the root such that the resulting matrix is traceless and has eigenvalues that obey Abs[lamda_i-lambda_j] <= 2 Pi.  The analogous Mathematica function MatrixLog[] uses a different convention for the branch taken in the complex plane.";
Options[SULog] = Options[getPhases];
SULog[mat_?MatrixQ, opts:OptionsPattern[]] :=
 Block[{phases, vectors, center},
  {phases, vectors, center} = getPhases[mat, True, opts];
  (* equivalent to Diagonal[I phases].Conjugate[vectors] *)
  Transpose[vectors].(
      (* equivalent to Diagonal[I phases].Conjugate[vectors] *)
      I phases * Conjugate[vectors])];
SUNorm::usage = "Distance of group element from the identity (or the nearest element of the center) using the tangent space. Returns the norm and the associated element of the center, labeled by associated element of Z_N.";
Options[SUNorm] = Options[getPhases];
SUNorm[mat_?MatrixQ, opts:OptionsPattern[]] :=
 Block[{phases, center},
   {phases, center} = getPhases[mat, False, opts];
   {Norm[phases], center}];


stringOperator::usage = "See arXiv:hep-lat/0107007v2 1 Aug 2001, Appendix A.  There is a notational error in Teper's Eqns. (49) and (50).";
stringOperator::unknown = "Unknown `1`";
stringOperator[uu_, {op1_, op2_}] :=
    {stringOperator[uu, op1], stringOperator[uu, op2]};
stringOperator[uu_, 1] := stringOperator[uu, "1"];
stringOperator[uu_, op_String]:=
 Block[{nc = Length[uu]},
   Which[
      (* Use the normalization from the Italian paper:
         action on the identity element is 1.
         In this case, the action on any group element lies inside
         the unit circle in the complex plane.
         For the antisymmetric representation dimension > nc,
         the action of the operator on any group element is zero. *)
      op == "1" || op == "t1", Tr[uu]/nc,
      op == "2S", (Tr[uu]^2 + Tr[uu.uu])/(nc (nc+1)),
      op == "2A", (Tr[uu]^2 - Tr[uu.uu])/(nc (nc-1)),
      op == "3S",
      (Tr[uu]^3 + 3 Tr[uu] Tr[uu.uu] + 2 Tr[uu.uu.uu])/(nc (nc+2) (nc+1)),
      op == "3A",
      If[nc > 2,
         (Tr[uu]^3 - 3 Tr[uu] Tr[uu.uu] + 2 Tr[uu.uu.uu])/(nc (nc-2) (nc-1)),
         0],
      op == "3M",
      (Tr[uu]^3 - Tr[uu.uu.uu])/(nc (nc+1) (nc-1)),
      op == "4S",
      (Tr[uu]^4 + 6 Tr[uu]^2 Tr[uu.uu] + 3 Tr[uu.uu]^2 +
       8 Tr[uu] Tr[uu.uu.uu] + 6 Tr[uu.uu.uu.uu])/(nc (nc+1) (nc+2) (nc+3)),
      op == "4A",
      If[nc > 3,
         (Tr[uu]^4 - 6 Tr[uu]^2 Tr[uu.uu] + 3 Tr[uu.uu]^2 +
          8 Tr[uu] Tr[uu.uu.uu] - 6 Tr[uu.uu.uu.uu])/
         (nc (nc-1) (nc-2) (nc-3)),
         0],
      (* various monomials *)
      op == "t2", Tr[uu.uu]/nc,
      op == "t1t1", Tr[uu]^2/nc^2,
      op == "t3", Tr[uu.uu.uu]/nc,
      op == "t2t1", Tr[uu.uu]*Tr[uu]/nc^2,
      op == "t1t1t1", Tr[uu]^3/nc^3,
      op == "t4", Tr[uu.uu.uu.uu]/nc,
      op == "t5", Tr[uu.uu.uu.uu.uu]/nc,
      op == "t6", Tr[uu.uu.uu.uu.uu.uu]/nc,
      op == "t7", Tr[uu.uu.uu.uu.uu.uu.uu]/nc,
      op == "t8", Tr[uu.uu.uu.uu.uu.uu.uu.uu]/nc,
      (* Complete set of physical observables.
        Returns a vector of nc reals.  All other quantities
        can be inferred from these quantities. *)
      op == "phases",
      Block[{phases = First[getPhases[uu, False]]},
            Sort[phases, Greater]],
      True,
      Message[stringOperator::unknown, op]; $Failed
  ]];


multiNorm::usage = "Norm squared of a rank-2 tensor in the adjoint \
representation projected onto various multiplets.  See Stefan Keppeler \
and Malin Sjödahl, \"Orthogonal multiplet bases in SU(N_c) color space,\"\
arXiv:1207.0609v2 [hep-ph] 2 Oct 2012.  This does not include the \
other 10 multiplet or separate out the additional multiplet that occurs \
for nc>3." 
multiNorm["singlet", mat_] := Tr[mat]; 
multiNorm["1", mat_] := Tr[mat]^2/(nc^2 - 1); 
multiNorm["8S", mat_] := 
 Block[{x = TensorContract[mat.SUSymmetric[], {{1, 2}}]}, 
  Total[x^2] nc/(nc^2 - 4)];
multiNorm["8A", mat_] := 
  Block[{x = TensorContract[mat.SUStructure[], {{1, 2}}]}, 
   Total[x^2]/nc];
multiNorm["27S", mat_] := ((Tr[mat.Transpose[mat]] + Tr[mat.mat])/2 - 
   multiNorm["8S", mat] - multiNorm[1, mat]); 
multiNorm["10", mat_] := ((Tr[mat.Transpose[mat]] - Tr[mat.mat])/2 - 
            multiNorm["8A", mat]);
multiNorm[j_?IntegerQ, mat_] := multiNorm[ToString[j], mat]; 


maxSURotation::usage = "Find a color rotation that maximizes the \
magnitude squared of the diagonal part of the logarithm of the given \
list of link fields."; 
maxSURotationF[uu_, vars_, weights_, nc_] := 
 Block[{vv = MatrixExp[I vars.SUGenerators[nc]]}, 
  Total[Map[
      Total[weights*Im[Diagonal[SULog[ConjugateTranspose[vv].#.vv]]]^2]&, 
           uu]]]/;VectorQ[vars, NumericQ]; 
maxSURotation[uu_] := 
  Block[{bb, nc = Length[First[uu]], result, vars, weights},
  weights = Table[1 - i/(100 nc), {i, nc}];
  vars = Array[bb, {nc^2 - 1}]; 
  result = FindMaximum[maxSURotationF[uu, vars, weights, nc], 
                       Map[{#, 0.0} &, vars]];
  If[False, Print["  result ", result]];
  MatrixExp[I (vars/.result[[2]]).SUGenerators[nc]]];
