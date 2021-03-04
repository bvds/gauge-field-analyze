(*
  Do one link at a time.

  That is, ignore off-diagonal color blocks of the Hessian.
  Numerically, this is much easier.
 *)

linkSaddlePointStep::usage = "One iteration of the saddle point search applied to a single link.  Set ignoreCutoff -> True to test this against the explicit calculation.  Returns the new link and the magnitude of the shift squared.";
Options[linkSaddlePointStep] = Join[
    {dampingFactor -> 1, printHessian -> False, printShift -> False,
    linkConcavity -> None},
    Options[applyCutoff1]];
linkSaddlePointStep[ff_?MatrixQ, opts:OptionsPattern[]] := 
 Block[{
   copts = Apply[Sequence, FilterRules[{opts}, Options[applyCutoff1]]],
   damping = OptionValue[dampingFactor],
   grad, mm, oo, values, deltas},
  grad = Map[Tr, SUGenerators[].ff];
  (* Solve linear system mm.shifts = -Im[grad],
    paying careful attention to nullspace of mm.
    Could equivalently use LDL^T factorization with symmetric pivoting. *)
  mm = Re[Tr[ff]]/(2.0 Length[ff]) IdentityMatrix[Length[grad]] + 
       SUSymmetric[].Re[grad]/2.0;
  {values, oo} = Eigensystem[mm];
  deltas = -damping applyCutoff1[values, oo.Im[grad], copts].oo;
  If[OptionValue[printHessian],
     Print["Hessian:  ", mm]; 
     Print["gradient:  ", Im[grad]]]; 
  If[OptionValue[printShift],
     Print["shifts:  ", deltas]];
  (* The 1/2 comes from our normalization convention for
    the generators of SU(N). *)
  {MatrixExp[I deltas.SUGenerators[]], Total[deltas^2]/2}
 ]/;OptionValue[linkConcavity] === None; 

linkSaddlePointStep[ff_?MatrixQ, opts:OptionsPattern[]] := 
 Block[{uu, damping = OptionValue[dampingFactor],
        sign = Sign[OptionValue[linkConcavity]],
        phases, vectors, center},
       uu = SUStapleMinimum[sign*ff];
       (* See the code for SUPower[] *)
       {phases, vectors, center} = getPhases[uu, True];
       {Transpose[vectors].(Exp[I phases damping] Conjugate[vectors]),
        Total[(phases damping)^2]}]/;NumberQ[OptionValue[linkConcavity]]; 

linkSaddlePoint::usage = "Saddle point search applied to a single link.  Returns new link value and norm of the shift squared.";
Options[linkSaddlePoint] = 
  Join[{Tolerance -> 10^-6, maxCount -> 1, printHessian -> False, 
        printDetails -> True},
       Options[linkSaddlePointStep]];
linkSaddlePoint[dir_, coords_, opts:OptionsPattern[]] := 
 If[boundaryLinkQ[dir, coords],
    {getLink[dir, coords], 0},
    Block[{uu = getLink[dir, coords], delta2 = 0,
           sopts = Apply[Sequence, FilterRules[
               {opts}, Options[linkSaddlePointStep]]],
           sumStaples = sumStaples[dir, coords],
           debug = printLevel[OptionValue[printDetails], 1]},
      Do[
          Block[{u1 = SUPower[uu, 0.5]},
                {uu, delta2} = linkSaddlePointStep[
                    u1.sumStaples.u1, sopts];
                uu = u1.uu.u1],
          OptionValue[maxCount]];
      If[debug > 0, SUMatrixQ[uu, Tolerance -> 10^-7]];
      {uu, delta2}]];

singleLinkStep::usage = "Apply one-link saddle point step to all links, updating in a checkerboard fashion, returning the size of the step.";
Options[singleLinkStep] = Options[linkSaddlePoint];
singleLinkStep[opts:OptionsPattern[]] :=
 Block[{norm2 = 0},
   Do[
       (* Shared variables are super-slow in Mathematica,
          so we find the new link values and then update
          gaugeField. *)
       Scan[
           (setLink[dir, #[[1]], #[[2, 1]]];
            norm2 += #[[2, 2]])&,
           ParallelTable[
               Block[{coords = latticeCoordinates[i]},
                 If[
                     Mod[Total[coords], 2] == checkerboard,
                     coords -> linkSaddlePoint[dir, coords, opts],
                     Nothing]],
               {i, latticeVolume[]}]
       ],
       {dir, nd}, {checkerboard, 0, 1}]; 
   Sqrt[norm2/(nd*latticeVolume[])]];


Options[makeObservableTrajectory] = Join[
    {"observables" -> {"distance", "norm", "polyakovCorrelator",
                       "wilsonDiagonal", "wilsonTriangle"},
     "polyakovLoopType" -> "simple",
     "skip" -> 1,
     "step" -> "../cache/3-3/step-16-28", 
     "periodic" -> "/mnt/samson-data/3-3/periodic-16-28",
     (* Defaults for the Newton's steps *)
     largeShiftOptions -> {eigenPairs -> -100, maxLanczosVecs -> 1000},
     Method -> {"External", maxLanczosVecs -> 1000},
     externalAction -> "detach", readInterval -> 5,
     remoteHost -> "samson", processes -> 8
    },
    Options[polyakovCorrelatorTallies],
    Options[singleLinkStep], Options[applyGaugeTransforms]]; 
makeObservableTrajectory[set_, label_, n_, 
                            opts:OptionsPattern[]] := 
 Block[{gaugeField, delta, distance = 0,
           shiftDistance,
           debug = printLevel[OptionValue[printDetails], 1],
        diagonalQ = StringMatchQ[label, "s*"],
        dd, results, gaugeSegments,
        gopts = Apply[Sequence, FilterRules[
            Join[{opts}, Options[makeObservableTrajectory]],
            Options[applyGaugeTransforms]]],
        popts = Apply[Sequence, FilterRules[
            {opts}, Options[polyakovCorrelatorTallies]]],
        sopts = Apply[Sequence, FilterRules[
            {opts}, Options[singleLinkStep]]]},
  Get[OptionValue["periodic"] <> "-" <> ToString[set] <> ".m"];
  
  results = Transpose[Table[
      Catch[
      If[debug > 0, Print["Starting step ", i]];
      gaugeSegments = Association[];
      If[i > 0,
         If[diagonalQ,
            dd = singleLinkStep[sopts],
            shiftDistance = Null; (* not always defined in the file *)
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
              observable == "averagePlaquette",
              (* This is redundant with wilsonDiagonal *)
              averagePlaquette[],
              observable == "norm",
              (* Except for distance, which is already accumulated above,
                the othere quantities calculated in this loop are gauge
                invariant, so fixing a gauge won't cause any harm. *)
              applyGaugeTransforms[Flatten[
                  {1, 1, 6, 1, Table[{5, 6}, {11}], 7, 6, 7, 6, 7}], gopts];
              latticeNorm[],
              observable == "polyakovCorrelator",
              Map[talliesToAverageErrors,
                  polyakovCorrelatorTallies[
                      OptionValue["polyakovLoopType"], "1", popts]],
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
              Print["Unknown observable ", observable]; Abort[]
              ]}, {observable, OptionValue["observables"]}], Nothing],
      Get::noopen],
      {i, 0, n}]];
  Block[{i=1}, Do[
      observableTrajectory[set, label, observable] = results[[i++]],
      {observable, OptionValue["observables"]}]]];

bulkPolyakovLoops::usage = "Aggregate Polyakov loop correlators over a number of lattice configurations.";
Options[bulkPolyakovLoops] = {"polyakovLoopType" -> "simple",
                              "lower" -> 1, "upper" -> 1};
bulkPolyakovLoops[inPath_, outFile_, opts:OptionsPattern[]] :=
 Block[{t0 = SessionTime[], t1, t2, t3, t4,
        tallyData, cov},
  tallyData = Merge[Table[
      Block[{gaugeField},
            Get[inPath <> ToString[i] <> ".m"]; 
            polyakovCorrelatorTallies[OptionValue["polyakovLoopType"], "1"]],
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
       stringModelBare[Map[talliesToAverageErrors, #], printResult -> False, 
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
    FilterRules[Options[stringModelBare], {Except["covarianceMatrix"]}]];
bulkStringModel[opts:OptionsPattern[]] :=
 Block[{sopts = Apply[Sequence, FilterRules[
        {opts}, Options[stringModelBare]]], result},
  result = Table[
      Block[{tallyData, ff},
            tallyData = Map[valueError[#[[i]], Null]&,
                                      polyakovCorrelatorValues];
            ff = stringModelBare[
                tallyData, sopts,
                "covarianceMatrix" -> polyakovCorrelatorCovariance];
            {i, ff["chiSquared"], Map[
                  {#[[1, 1]], valueError[#[[1, 2]], #[[2]]]}&, 
                     Transpose[
                         ff[{"BestFitParameters", 
                             "ParameterErrors"}]]]}],
      {i, OptionValue["lower"], OptionValue["upper"]}];
  Print["Delete graphic before saving notebook."];
  Print[Graphics[
      Map[{Text[#[[1]], {#[[2]], #[[3, 1, 2, 1]]}], 
           Line[{{#[[2]], #[[3, 1, 2, 1]] - #[[3, 1, 2, 2]]},
                 {#[[2]], #[[3, 1, 2, 1]] + #[[3, 1, 2, 2]]}}]}&,
      result], Frame -> True, 
      AspectRatio -> 1, PlotRange -> {All, {0, 0.030}}]];
  result];

bulkWilsonLoops::usage = "Aggregate Wilson loops over a number of lattice configurations.  Skip some sizes, to save on time.";
Options[bulkWilsonLoops] = {"lower" -> 1, "upper" -> 1,
                            "lowerCutoff" -> 1, "skip" -> 1};
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
