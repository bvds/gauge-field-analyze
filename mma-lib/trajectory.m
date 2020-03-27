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
  u1.MatrixExp[I deltas.oo.SUGenerators[]]]; 

linkSaddlePoint::usage = "Saddle point search applied to a single link.  Returns new link value and norm of the change.";
Options[linkSaddlePoint] = 
  Join[{Tolerance -> 10^-6, maxCount -> 1, printHessian -> False, 
        printDetails -> True},
       Options[linkSaddlePointStep]];
linkSaddlePoint[dir_, coords_, opts:OptionsPattern[]] := 
 Block[{maxShift = 4 Pi, u0 = getLink[dir, coords], u1, u2,
        linkShiftNorm,
        sopts = Apply[Sequence, FilterRules[
            {opts}, Options[linkSaddlePointStep]]],
        u2SumStaples, count = 0, uu,
        debug = printLevel[OptionValue[printDetails], 1]},
  u2 = u1 = SUPower[u0, 0.5];
  u2SumStaples = u2.sumStaples[dir, coords];
  While[maxShift > OptionValue[Tolerance] && 
	count < OptionValue[maxCount],
	count++;
	u1 = linkSaddlePointStep[u2SumStaples, u1, sopts]];
  uu = u1.u2;
  If[debug > 0, SUMatrixQ[uu, Tolerance -> 10^-7]];
  {uu, First[SUNorm[uu.ConjugateTranspose[u0]]]}];

makeCoordList[] := Block[
    {result = Array[Null&, latticeVolume[]]},
    Do[Block[{coord = latticeCoordinates[k]},
             result[[linearSiteIndex[coord]]] = coord],
       {k, latticeVolume[]}]; result];
singleLinkStep::usage = "Apply one-link saddle point step to all links, returning the updated links and the size of the step.";
Options[singleLinkStep] = Options[linkSaddlePoint];
singleLinkStep[coordList_List:Null, opts:OptionsPattern[]] :=
    Block[{tt = Table[
              ParallelMap[
                  linkSaddlePoint[dir, #, opts]&,
                                 If[coordList === Null,
                                    makeCoordList[], coordList]],
	      {dir, nd}]},
          {Map[First, tt, {2}],
           Sqrt[Mean[Flatten[Map[Last[#]^2&, tt, {2}]]]]}];


Options[makeObservableTrajectory] = Join[
    Options[singleLinkStep], {
        "observables" -> {"distance", "polyakovCorrelator",
                          "wilsonDiagonal", "wilsonTriangle"},
        "skip" -> 1,
        "step" -> "../data/3-3/step-16-28", 
        "periodic" -> "../data/3-3/periodic-16-28"}]; 
makeObservableTrajectory[set_, label_, n_, 
                            opts:OptionsPattern[]] := 
 Block[{gaugeField, delta, gaugeField0, distance = 0,
        lastGaugeField, coordList, dd, results,
        sopts = Apply[Sequence, FilterRules[
            {opts}, Options[singleLinkStep]]]},
  Get[OptionValue["periodic"] <> "-" <> ToString[set] <> ".m"];
  If[StringMatchQ[label, "single*"],
     coordList = makeCoordList[]];
  lastGaugeField = gaugeField;
  
  results = Transpose[Table[
      Print["Starting ", i]; 
      If[i > 0,
         If[StringMatchQ[label, "single*"],
            {gaugeField, dd} = singleLinkStep[coordList, sopts],
            Get[StringRiffle[{OptionValue["step"], ToString[set],
                              label, ToString[i]},"-"]<>".m"];
            gaugeField = gaugeField0;
            dd = If[True,
                    (* These are equivalent *)
                    Norm[delta]/Sqrt[2*nd*latticeVolume[]],
                    latticeDistance[gaugeField, lastGaugeField];
                    lastGaugeField = gaugeField]];
         distance += dd];
      If[Mod[i, OptionValue["skip"]] == 0, Table[
          {i, Which[
              observable == "distance",
              distance,
              observable == "polyakovCorrelator",
              talliesToAverageErrors[
                  polyakovCorrelatorTallies["simple","1"]],
              observable == "wilsonDiagonal",
              Merge[Apply[averageWilsonLoop,
                           Table[{w, w}, {w, Max[latticeDimensions]-1}],
                           {1}], Total],
              observable == "wilsonTriangle",
              Merge[Apply[averageWilsonLoop,
                           Flatten[Table[{w1, w2},
                                         {w1, 4, Max[latticeDimensions]-1,  2},
                                         {w2, 4, w1, 2}], 1], {1}], Total],
              True,
              Abort[]
              ]}, {observable, OptionValue["observables"]}], Nothing],
      {i, 0, n}]];
  Block[{i=1}, Do[
      observableTrajectory[set, label, observable] = results[[i++]],
      {observable, OptionValue["observables"]}]]];

bulkPolyakovLoops::usage = "Aggregate Polyakov loop correlators over a number of lattice configurations.";
bulkPolyakovLoops[inPath_, outFile_, steps_] :=
Block[{tallyData = Table[Block[{gaugeField}, 
        Get[inPath <> ToString[i] <> ".m"]; 
        polyakovCorrelatorTallies["simple", "1"]], {i, 1, steps}],
       data}, 
 configurationCount = Length[tallyData];
 data = Transpose[Map[#[[2]]/#[[1]] &, Values[tallyData], {2}]];
 covarianceEigenvalues = Table[
     Block[{dd = Map[Take[#, k]&, data]},
           Map[{k, #}&, Chop[Sort[Eigenvalues[
               Outer[Covariance, dd, dd, 1]]]]]],
     {k, 2, configurationCount}];
 polyakovCorrelatorCovariance = Outer[Covariance, data, data, 1];
 stringModelValuesState = 
  Map[Block[{ff = 
       stringModel[talliesToAverageErrors[#], printResult -> False, 
        "eigenstate" -> 2, "constantTerm" -> False]}, 
     Append[Map[valueError[#[[1]], #[[2]]] &, 
       ff["ParameterTableEntries"]], {ff["chiSquared"], 
       Length[#] - Length[ff["BestFitParameters"]]}]] &, tallyData]; 
 stringModelValuesPotential = 
  Map[Block[{ff = 
       stringModel[talliesToAverageErrors[#], printResult -> False, 
        "eigenstate" -> 0, "order" -> 0, "constantTerm" -> False]}, 
     Append[Map[valueError[#[[1]], #[[2]]] &, 
       ff["ParameterTableEntries"]], {ff["chiSquared"], 
       Length[#] - Length[ff["BestFitParameters"]]}]] &, tallyData]; 
 polyakovCorrelatorMerged = 
  talliesToAverageErrors[Merge[tallyData, Total]];
 polyakovCorrelatorGrouped = 
  talliesToAverageErrors[
   Merge[Map[{1, #[[2]]/#[[1]], (#[[2]]/#[[1]])^2} &, tallyData, {2}],
         Total]];
 DeleteFile[outFile];
 Save[outFile,
      {"polyakovCorrelatorMerged", "polyakovCorrelatorGrouped",
       "polyakovCorrelatorCovariance", "configurationCount",
       "covarianceEigenvalues",
       "stringModelValuesState", "stringModelValuesPotential"}]]
