Options[makeStringModelTrajectory] = {
    "step" -> "../data/3-3/step-16-28-5", 
    "periodic" -> "../data/3-3/periodic-16-28-5", "lowerCutoff" -> 40}; 
makeStringModelTrajectory[label_, n_, OptionsPattern[]] := 
 Block[{gaugeField, tallyData, ff, gaugeField0, distance = 0, 
   lastGaugeField}, Get[OptionValue["periodic"] <> ".m"]; 
  lastGaugeField = gaugeField;
  stringModelTrajectory[label] = 
   Table[Print["Starting ", i]; 
    If[i > 0, 
     Get[OptionValue["step"] <> "-" <> label <> "-" <> ToString[i] <> 
       ".m"]; gaugeField = gaugeField0;
     distance += latticeDistance[gaugeField, lastGaugeField]; 
     lastGaugeField = gaugeField]; 
    tallyData = talliesToAverageErrors[polyakovLoopTallies["simple"]];
     ff = stringModel[tallyData, 
      "lowerCutoff" -> OptionValue["lowerCutoff"]]; 
    Join[{i, distance}, 
     MapThread[
      valueError[#1[[2]], #2] &, {ff["BestFitParameters"], 
                                  ff["ParameterErrors"]}]], {i, 0, n}]];

Options[makeWilsonTrajectory] = {
    "step" -> "../data/3-3/step-16-28-5",
    "periodic" -> "../data/3-3/periodic-16-28-5"}; 
makeWilsonTrajectory[label_, n_, OptionsPattern[]] := 
 Block[{gaugeField, gaugeField0, distance, lastGaugeField}, 
  Get[OptionValue["periodic"] <> ".m"];
  Print["String tension ", 
   stringTension = teperTension[3, 3, 1, beta, averagePlaquette[]]^2];
   wilsonTrajectory[label] = 
   Table[Print["Starting ", i]; 
    If[i > 0, 
     Get[OptionValue["step"] <> "-" <> label <> "-" <> ToString[i] <> 
       ".m"]; gaugeField = gaugeField0;
     distance += latticeDistance[gaugeField, lastGaugeField]; 
     lastGaugeField = gaugeField]; {i, distance, 
     averageWilsonLoop[1, 1], averageWilsonLoop[2, 2], 
     averageWilsonLoop[4, 4], averageWilsonLoop[6, 6], 
     averageWilsonLoop[8, 8]}, {i, 0, n}]];

(*
  Do one link at a time.

  That is, ignore off-diagonal color blocks of the Hessian.
  Numerically, this is much easier.
 *)

linkSaddlePointStep::usage = "One iteration of the saddle point search applied to a single link.  Set ignoreCutoff -> True to test this against the explicit calculation.";
Options[linkSaddlePointStep] = {ignoreCutoff -> False,
    linkCutoff -> 1, rescaleCutoff -> 1, dampingFactor -> 1,
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
makeCoordList[] := Block[
    {result = Array[Null&, latticeVolume[]]},
    Do[Block[{coord = latticeCoordinates[k]},
             result[[linearSiteIndex[coord]]] = coord],
       {k, latticeVolume[]}]; result];
makeSingleLinkSaddlepointTrajectory::usage = "Apply one-link saddle point step to all links, returning associated statistics.  This starts with some given gaugeField and updates it.";
Options[makeSingleLinkSaddlepointTrajectory] = Join[
    {maxAvgPlaquette -> 1, "periodic" -> "../data/3-3/periodic-16-28-5"},
    Options[linkSaddlePoint]];
makeSingleLinkSaddlepointTrajectory[n_, opts:OptionsPattern[]] :=
 Block[{tallyData, avgP = 0, ff, result, t0, t1,
        coordList = makeCoordList[],
        lspo = Apply[Sequence, FilterRules[{opts}, Options[linkSaddlePoint]]],
        debug = printLevel[OptionValue[printDetails], 1]},
  gfi = gaugeField;
  (* Sow[] and Reap[] allow an early exit from the loop. *)
  Reap[Do[
      If[avgP >= OptionValue[maxAvgPlaquette], Break[]];
      t0 = SessionTime[];
      If[i>0,
         gaugeField = Table[
             ParallelMap[First[linkSaddlePoint[dir, #, lspo]]&,
                              coordList],
	     {dir, nd}];
         ];
      t1 = SessionTime[];
      If[debug > 1, Print["Step ", i, " in ", t1 - t0, " seconds"]];
      tallyData = talliesToAverageErrors[polyakovLoopTallies["smeared"]];
      ff = exponentialModel[tallyData];
      avgP = averagePlaquette[];
      result = Join[{i, avgP},
           Inner[{#1[[1]], #1[[2]], #2} &, ff["BestFitParameters"], 
                 ff["ParameterErrors"], List]];
      If[debug > 0, Print[result]];
      Sow[result],
      {i, 0, n}]][[2,1]]];