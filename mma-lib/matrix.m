(*
  Utilities for sparse matrix construction;
  used by saddle-point.m and gauge.m.
 *)
printMemory[n_] := EngineeringForm[N[ByteCount[n]], 3];
printLevel[opt_, default_] :=
    Which[NumberQ[opt], opt,
          opt === True, default, opt === False, 0,
          True, default];

symAdd::usage = oneAdd::useage = "If full=False, take the lower triangle.  The export to an *.mtx file dumps the lower triangle of a symmetric matrix.";
symAdd::index = "Invalid index.";
SetAttributes[symAdd, HoldAll];
symAdd[m_, i_, j_, value_, full_:True] :=
    (If[full || i>=j, m[{i, j}] = Lookup[m, Key[{i, j}], 0.0] + value];
     If[full || j>i, m[{j, i}] = Lookup[m, Key[{j, i}], 0.0] + value]);
SetAttributes[oneAdd, HoldAll];
oneAdd[m_, i_, j_, value_, full_:True] :=
    If[full || i>=j, m[{i, j}] = Lookup[m, Key[{i, j}], 0.0] + value];
reTrDot::usage = "Compute Re[Tr[a.b]] where a and b are color matrices.";
reTrDot = Compile[{{x, _Complex, 2}, {y, _Complex, 2}},
	Sum[Re[x[[i, j]] y[[j, i]]], {i, nc}, {j, nc}], {{nc, _Integer}}];
