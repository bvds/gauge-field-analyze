(*
  Routines to expedite timing of a block of code
 *)
SetAttributes[addTime, {HoldAll, SequenceHold}];
addTime[timer_, expr_] :=
 Block[{t1 = SessionTime[], result = expr},
       If[NumberQ[timer], timer += SessionTime[] - t1]; result];

SetAttributes[addTimeNull, {HoldAll, SequenceHold}];
addTimeNull[timer_, expr_] :=
    Block[{t1 = SessionTime[]},
          expr; If[NumberQ[timer], timer += SessionTime[] - t1]];
