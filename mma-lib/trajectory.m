(*
  Do one link at a time.

  That is, ignore off-diagonal color blocks of the Hessian.
  Numerically, this is much easier.
 *)

linkSaddlePointStep::usage = "One iteration of the saddle point search applied to a single link.  Set ignoreCutoff -> True to test this against the explicit calculation.";
Options[linkSaddlePointStep] = {ignoreCutoff -> False,
    linkCutoff -> 1, dampingFactor -> 1,
    printHessian -> False, printShift -> False};
linkSaddlePointStep[u2Staples_, u1_, OptionsPattern[]] := 
 Block[{
     cutoff = If[OptionValue[ignoreCutoff], 1, OptionValue[linkCutoff]],
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
 If[boundaryLinkQ[dir, coords],
    {getLink[dir, coords], 0},
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
     {uu, First[SUNorm[uu.ConjugateTranspose[u0]]]}]];

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
    {"observables" -> {"distance", "norm", "polyakovCorrelator",
                       "wilsonDiagonal", "wilsonTriangle"},
     "skip" -> 1,
     "step" -> "../cache/3-3/step-16-28", 
     "periodic" -> "/mnt/samson-data/3-3/periodic-16-28",
     (* Defaults for the Newton's steps *)
     largeShiftOptions -> {eigenPairs -> -100, maxLanczosVecs -> 1000},
     Method -> "External", externalAction -> "detach", readInterval -> 5,
     remoteHost -> "samson", processes -> 8
    },
    Options[singleLinkStep], Options[applyGaugeTransforms]]; 
makeObservableTrajectory[set_, label_, n_, 
                            opts:OptionsPattern[]] := 
 Block[{gaugeField, delta, distance = 0,
           shiftDistance,
           debug = printLevel[OptionValue[printDetails], 1],
        diagonalQ = StringMatchQ[label, "s*"],
        coordList, dd, results, gaugeSegments,
        gopts = Apply[Sequence, FilterRules[
            Join[{opts}, Options[makeObservableTrajectory]],
            Options[applyGaugeTransforms]]],
        sopts = Apply[Sequence, FilterRules[
            {opts}, Options[singleLinkStep]]]},
  Get[OptionValue["periodic"] <> "-" <> ToString[set] <> ".m"];
  If[diagonalQ,
     coordList = makeCoordList[]];
  
  results = Transpose[Table[
      Catch[
      If[debug > 0, Print["Starting step ", i]];
      gaugeSegments = Association[];
      If[i > 0,
         If[diagonalQ,
            {gaugeField, dd} = singleLinkStep[coordList, sopts],
            shiftDistance = Null; (* not always defined *)
            Check[
                Get[StringRiffle[{OptionValue["step"], ToString[set],
                                  label, ToString[i]}, "-"]<>".m"],
                Print["Skipping ", i, ":  file missing."];
                Throw[Nothing, Get::noopen],
                {Get::noopen}];
            (* Delta is not defined for the
              diagonal-only steps. *)
            dd = If[
                shiftDistance =!= Null,
                shiftDistance,
                Norm[delta]/Sqrt[2*nd*latticeVolume[]]]];
         distance += dd];
      If[Mod[i, OptionValue["skip"]] == 0, Table[
          {i, If[debug > 0, Print["  ", DateObject[], " start ", observable]];
              Which[
              observable == "distance",
              distance,
              observable == "norm",
              (* Except for distance, which is already accumulated above,
                the quantities are gauge invariant, so fixing a gauge
                won't cause any harm. *)
              applyGaugeTransforms[Flatten[
                  {1, 1, 6, 1, Table[{5, 6}, {12}], 7, 7, 7}], gopts];
              latticeNorm[],
              observable == "polyakovCorrelator",
              Map[talliesToAverageErrors,
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
      Get::noopen],
      {i, 0, n}]];
  Block[{i=1}, Do[
      observableTrajectory[set, label, observable] = results[[i++]],
      {observable, OptionValue["observables"]}]]];

bulkPolyakovLoops::usage = "Aggregate Polyakov loop correlators over a number of lattice configurations.";
Options[bulkPolyakovLoops] = {"lower" -> 1, "upper" -> 1};
bulkPolyakovLoops[inPath_, outFile_, opts:OptionsPattern[]] :=
 Block[{t0 = SessionTime[], t1, t2, t3, t4,
        tallyData, cov},
  tallyData = Merge[Table[
      Block[{gaugeField},
            Get[inPath <> ToString[i] <> ".m"]; 
            polyakovCorrelatorTallies["simple", "1"]],
      {i, OptionValue["lower"], OptionValue["upper"]}], Join];
  t1 = SessionTime[];
  Print["Finished aggregating data. tallyData memory:  ",
        ByteCount[tallyData], " in ", t1-t0, " s"];
  polyakovCorrelatorValues = Map[#[[2]]/#[[1]] &, tallyData, {2}];
  (* Only include pairs of Polyakov loop correlators that
    can be in the same direction. See analysis in gauge.nb *)
  polyakovCorrelatorCovariance =
   Block[{data = Normal[polyakovCorrelatorValues]},
     Outer[If[#1[[1, -1]] == #2[[1, -1]], Covariance[#1[[2]], #2[[2]]], 0]&,
             data, data, 1]];
  (* version for one-configuration fit. *)
  cov = If[Min[Eigenvalues[polyakovCorrelatorCovariance]] > 0,
           polyakovCorrelatorCovariance,
           (* One could revert back to "None" but this is better
             for code testing before scaling up to a longer
             calculation. *)
           Message[General::npdef, polyakovCorrelatorCovariance];
           DiagonalMatrix[Diagonal[polyakovCorrelatorCovariance]]];
  t2 = SessionTime[];
  Print["Finished correlators in ", t2-t1, "s"];
  stringModelValuesState = Null;
  stringModelValuesPotential = Null;
  (* stringModelValuesState = 
   Map[Block[{ff = 
       stringModel[Map[talliesToAverageErrors, #], printResult -> False, 
                   "covarianceMatrix" -> cov,
                   "eigenstate" -> 2]},
     Append[Map[valueError[#[[1]], #[[2]]] &, 
       ff["ParameterTableEntries"]], {ff["chiSquared"], 
       Length[#] - Length[ff["BestFitParameters"]]}]] &, tallyData]; 
 stringModelValuesPotential = 
  Map[Block[{ff = 
       stringModel[Map[talliesToAverageErrors, #], printResult -> False, 
                   "covarianceMatrix" -> cov,
                   "eigenstate" -> 0, "order" -> 0]}, 
     Append[Map[valueError[#[[1]], #[[2]]] &, 
       ff["ParameterTableEntries"]], {ff["chiSquared"], 
       Length[#] - Length[ff["BestFitParameters"]]}]] &, tallyData]; *)
  t3 = SessionTime[];
  Print["Finished models in ", t3-t2, " s"];
  polyakovCorrelatorMerged = Map[
      talliesToAverageErrors[Total[#]]&, tallyData];
  polyakovCorrelatorGrouped = Map[
  talliesToAverageErrors[
      Total[Map[{1, #, #^2} &, #]]]&,
                              polyakovCorrelatorValues];
  DeleteFile[outFile];
  Save[outFile,
       {"polyakovCorrelatorMerged", "polyakovCorrelatorGrouped",
        "polyakovCorrelatorCovariance", "polyakovCorrelatorValues",
        "stringModelValuesState", "stringModelValuesPotential"}];
  t4 = SessionTime[];
  Print["Merged and grouped correlators, dump results ", t4-t3, " s"]];

bulkStringModel::usage = "Aggregate string model fits over a number of lattice configurations.";
Options[bulkStringModel] = Join[
    {"lower" -> 1, "upper" -> 1},
    FilterRules[Options[stringModel], {Except["covarianceMatrix"]}]];
bulkStringModel[opts:OptionsPattern[]] :=
 Block[{sopts = Apply[Sequence, FilterRules[
        {opts}, Options[stringModelFit]]]},
  Table[
      Block[{tallyData, ff},
            tallyData = Map[valueError[#[[i]], Null]&,
                                      polyakovCorrelatorValues];
            ff = stringModel[
                tallyData, sopts,
                "covarianceMatrix" -> polyakovCorrelatorCovariance];
            {i, ff["chiSquared"], Map[
                  {#[[1, 1]], valueError[#[[1, 2]], #[[2]]]}&, 
                     Transpose[
                         ff[{"BestFitParameters", 
                             "ParameterErrors"}]]]}],
      {i, OptionValue["lower"], OptionValue["upper"]}]];

bulkWilsonLoops::usage = "Aggregate Wilson loops over a number of lattice configurations.  Skip some sizes, to save on time.";
Options[bulkWilsonLoops] = {"lower" -> 1, "upper" -> 1,
                            "lowerCutoff" -> 4, "skip" -> 2};
bulkWilsonLoops[inPath_, outFile_, opts:OptionsPattern[]] :=
 Block[{t0 = SessionTime[], t1, t2, t3,
        tallyData, cov},
  (* Parallelization must be implemented at the gaugeField level so
    that the gaugeSegments cache works properly. *)
  tallyData = Merge[ParallelTable[
      Block[{gaugeField, gaugeSegments = Association[]},
            Get[inPath <> ToString[i] <> ".m"]; 
            Merge[Flatten[
                Table[wilsonLoopTallies[w1, w2],
                      {w1, OptionValue["lowerCutoff"],
                       Max[latticeDimensions] - 1, OptionValue["skip"]},
                      {w2, OptionValue["lowerCutoff"], w1,
                       OptionValue["skip"]}]],
                  Total]],
      {i, OptionValue["lower"], OptionValue["upper"]}], Join];
  t1 = SessionTime[];
  Print["Finished aggregating data. tallyData memory:  ",
        ByteCount[tallyData], " in ", t1-t0, " s"];
  wilsonTriangleValues = Map[#[[2]]/#[[1]] &, tallyData, {2}];
  (* Only include pairs of Polyakov loop correlators that
    can be in the same direction. See analysis in gauge.nb *)
  wilsonTriangleCovariance =
  Block[{data = Normal[wilsonTriangleValues]},
        Outer[If[Drop[First[#1], 2] == Drop[First[#2], 2],
                 Covariance[#1[[2]], #2[[2]]], 0]&,
                data, data, 1]];
  (* version for one-configuration fit. *)
  cov = If[Min[Eigenvalues[wilsonTriangleCovariance]] > 0,
           wilsonTriangleCovariance,
           (* One could revert back to "None" but this is better
             for code testing before scaling up to a longer
             calculation. *)
           Message[General::npdef, wilsonTriangleCovariance];
           DiagonalMatrix[Diagonal[wilsonTriangleCovariance]]];
  t2 = SessionTime[];
  Print["Finished correlators in ", t2-t1, "s"];
  wilsonTriangleMerged = Map[
      talliesToAverageErrors[Total[#]]&, tallyData];
  wilsonTriangleGrouped = Map[
      talliesToAverageErrors[
          Total[Map[{1, #, #^2} &, #]]]&,
                            wilsonTriangleValues];
  DeleteFile[outFile];
  Save[outFile,
       {"wilsonTriangleMerged", "wilsonTriangleGrouped",
        "wilsonTriangleCovariance", "wilsonTriangleValues"}];
  t3 = SessionTime[];
  Print["Merged and grouped correlators, dump results ", t3-t2, " s"]];
