(*
  Do one link at a time.

  In this case, the resulting action density is dominated by
  one-plaquette spikes.
 *)

linkSaddlePointStep::usage = "One iteration of the saddle point search applied to a single link.  Set ignoreCutoff -> True to test this against the explicit calculation.";
Options[linkSaddlePointStep] = {ignoreCutoff -> False,
    linkCutoff :> Sqrt[nc^2/(3 beta)], rescaleCutoff -> 1, dampingFactor -> 1,
    printHessian -> False, printShift -> False};
linkSaddlePointStep[u2Staples_, u1_, OptionsPattern[]] := 
 Block[{
     cutoff = If[OptionValue[ignoreCutoff], 1,
                 OptionValue[linkCutoff]*OptionValue[rescaleCutoff]],
     zzz = If[OptionValue[ignoreCutoff], Infinity, 1],
   damping = OptionValue[dampingFactor], debug = False, 
   ff = u2Staples.u1, grad, mm, oo, values, deltas},
  grad = Map[Tr, SUGenerators[].ff];
  (* Solve linear system mm.shifts = -Im[grad],
    paying careful attention to nullspace of mm.
    Could equivalently use LDL^T factorization with symmetric pivoting. *)
  mm = Re[Tr[ff]]/(2.0 Length[ff]) IdentityMatrix[Length[grad]] + 
       SUSymmetric[].Re[grad]/2.0;
  {values, oo} = Eigensystem[mm];
  deltas = -damping applyCutoff1[values, oo.Im[grad], cutoff, zzz];
  If[OptionValue[printHessian],
     Print["Hessian:  ", mm]; 
     Print["gradient:  ", Im[grad]]]; 
  If[OptionValue[printShift],
     Print["shifts:  ", deltas.oo]];
  {u1.MatrixExp[I deltas.oo.SUGenerators[]], Max[Abs[deltas]]}]; 

linkSaddlePoint::usage = "Saddle point search applied to a single link.";
Options[linkSaddlePoint] = 
  Join[{Tolerance -> 10^-6, maxCount -> 1, printHessian -> False, 
        printDetails -> True},
       Options[linkSaddlePointStep]];
linkSaddlePoint[dir_, coords_, opts : OptionsPattern[]] := 
 Block[{maxShift = 4 Pi, u1 = SUPower[getLink[dir, coords], 0.5], u2, 
        u2SumStaples, count = 0, uu,
        debug = printLevel[OptionValue[printDetails], 1]},
  u2 = u1;
  u2SumStaples = u2.sumStaples[dir, coords];
  While[maxShift > OptionValue[Tolerance] && 
	count < OptionValue[maxCount],
	count++;
	{u1, maxShift} = linkSaddlePointStep[
	    u2SumStaples, u1, Apply[Sequence, 
	     FilterRules[{opts}, Options[linkSaddlePointStep]]]]; 
	If[debug > 2, Print["Max shift:  ", maxShift]]]; 
  uu = u1.u2;
  If[debug > 0, SUMatrixQ[uu, Tolerance -> 10^-7]];
  {uu, maxShift}];
allLinkSaddlePoint::usage = "Apply one-link saddle point step to all links, returning associated statistics.  This starts with some given gaugeField and updates it.";
Options[allLinkSaddlePoint] = Options[linkSaddlePoint];
allLinkSaddlePoint[n_, opts:OptionsPattern[]] :=
 Block[{distance, tallyData, ff, result, t0, t1,
        coordList = Sort[Table[latticeCoordinates[k], {k, latticeVolume[]}],
                         Order[linearSiteIndex[#1], linearSiteIndex[#2]]&],
        debug = printLevel[OptionValue[printDetails], 1]},
  distance = latticeDistanceFrom[gaugeField];
  Table[
      t0 = SessionTime[];
      If[i>0,
         gaugeField = Table[
             ParallelMap[First[linkSaddlePoint[dir, #, opts]]&,
                              coordList],
	     {dir, nd}];
         ];
      t1 = SessionTime[];
      If[debug > 1, Print["Step ", i, " in ", t1 - t0, " seconds"]];
      tallyData = talliesToAverageErrors[polyakovLoopTallies["smeared"]];
      ff = exponentialModel[tallyData];
      result = Join[{i, distance[gaugeField], averagePlaquette[]},
           Inner[{#1[[1]], #1[[2]], #2} &, ff["BestFitParameters"], 
                 ff["ParameterErrors"], List]];
      If[debug > 0, Print[result]];
      result,
      {i, 0, n}]];
