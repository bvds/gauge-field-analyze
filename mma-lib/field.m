opsList2 = {"2l", "2t1", "2t2"};
opsList3 = {"3pt1", "3pt2", "3m1a", "3m1b", "3m2", 
            "3m3", "3m4", "3m5", "3m6", "3m7", "3m8", "3m9", "3m10"};
componentMultiplet3[dir1_, dir2_, d1_, d3_, d2_] :=
 Block[{
   mid3 = Association[{
        {"long", "long", "long"} -> {"3m1b", Null},
        {"long", "long", "trans"} -> {"3m2", Null},
        {"long", "trans", "long"} -> {"3m3", Null},
        {"long", "trans", "trans"} -> {"3m4", Null},
        {"long", "out", "out"} -> {"3m5", Null},
        {"trans", "long", "long"} -> {"3m2", {1, -1, -1}},
        {"trans", "long", "trans"} -> {"3m6", Null},
        {"trans", "trans", "long"} -> {"3m4", {1, 1, -1}},
        {"trans", "trans", "trans"} -> {"3m7", Null},
        {"trans", "out", "out"} -> {"3m8", Null},
        {"out", "long", "out"} -> {"3m9", Null},
        {"out", "trans", "out"} -> {"3m10", Null},
        {"out", "out", "long"} -> {"3m5", {1, 1, -1}},
        {"out", "out", "trans"} -> {"3m8", {1, 1, 1}}}],
   key = Map[
       Which[# == dir1, "long", # == dir2, "trans", True, "out"]&,
       {d1, d3, d2}]},
       Lookup[mid3, Key[key], {Null, Null}]]/;dir1!=dir2;

linkCorrPlot3::usage = "Helper for plot routines, returning:  aspect ratio, sign of the data, and axis labels.";
linkCorrPlot3[] :=
 Block[{label1 = "#\[CapitalDelta]",
        labelp = "#\[CapitalDelta]", label2 = "r",
        label3 = "cov", labelm = "-cov"},
  Association[{
      "3pt1" -> {1, 1, {label1, labelp, label3}},
      "3pt2" -> {1/2, -1, {label1, labelp, label4}},
      "3m1a" -> {Sqrt[1 + (#[[-2]]/#[[-1]])^2]&[latticeDimensions],
                 1, {label1, label2, label3}},
      "3m1b" -> {1, 1, {label1, label2, label3}},
      "3m2" -> {1, 1, {label1, label2, label3}}, 
      "3m3" -> {1, -1, {label1, label2, labelm}},
      "3m4" -> {1, 1, {label1, label2, label3}},
      "3m5" -> {1, 1, {label1, label2, label3}},
      "3m6" -> {1, -1, {label1, label2, labelm}},
      "3m7" -> {1, 1, {label1, label2, label3}},
      "3m8" -> {1, -1, {label1, label2, labelm}},
      "3m9" -> {1, -1, {label1, label2, labelm}},
      "3m10" -> {1, 1, {label1, label2, label3}}}]];


Options[linkCorrelators] = {"symmetric" -> False,
    "distFlag" -> True, "operators" -> All, "timing" -> False};
linkCorrelators[n_, OptionsPattern[]] :=
 Block[{
    distFlag = OptionValue["distFlag"],
    ops = If[OptionValue["operators"] === All,
             If[n==2, opsList2, opsList3], OptionValue["operators"]],
    alf = Map[2*Im[Map[Tr, SUGenerators[].SULog[#]]] &, gaugeField, {2}], 
    t0 = SessionTime[], t1, result,
    half = IdentityMatrix[nd]/2, 
    vector = Function[{dir, len}, Table[If[i == dir, len, 0], {i, nd}]],
    cross = Function[{k1, k2}, If[k1 != k2, 6 - k1 - k2, $Failed]],
    displace = Function[dx, 
       Mod[dx + latticeDimensions/2, latticeDimensions] -
       latticeDimensions/2],
    fabc = N[If[OptionValue["symmetric"], SUSymmetric[], SUStructure[]]],
    tensorSign = If[OptionValue["symmetric"], 1, -1],
    addTimeNull = If[
        OptionValue["timing"],
        Function[{timer, expr},
                 Block[{t1 = SessionTime[]},
                       expr; timer += SessionTime[] - t1],
                 {HoldAll, SequenceHold}],
        Null&],
    update2, update3},
   update2 = Function[
       {tallies, op, key, w1, w2}, 
       If[KeyExistsQ[tallies, op], 
          Block[{lt = Lookup[tallies[op], Key[key],
                             If[distFlag,
                                {},
                                Join[{0, 0.0, 0.0},
                                     Array[0.0&, {2, Length[w1]}]]]],
                w12 = w1.w2},
                If[distFlag,
                   lt = {lt, w12},
                   lt += {1, w12, w12^2, w1, w2}];
                tallies[op][key] = lt]], {HoldFirst}];
   update3 = Function[
       {tallies, op, key, w1, w2, w3}, 
       If[KeyExistsQ[tallies, op], 
          Block[{f321, f32, f1, lt = Lookup[
              tallies[op], Key[key],
              If[distFlag, {},
                 Join[{0, 0.0, 0.0}, Array[0.0&, {6, Length[w1]}]]]]},
                If[distFlag,
                   lt = {lt, ((fabc.w3).w2).w1},
                   (* Right multiplication is faster for SparseArray.
                     Use the cyclic property of dabc or fabc. *)
                   lt += {1, f321=(f32=(fabc.w3).w2).w1, f321^2,
                          w1, w2, w3,
                          f32, (f1=fabc.w1).w3, tensorSign*f1.w2}];
                tallies[op][key] = lt]], {HoldFirst}];
   Print["ops: ", ops];
   If[False, Print[TableForm[
       Transpose[
           Map[Map[valueError[Mean[#], 
            StandardDeviation[#]/Sqrt[Length[#]]] &, Transpose[#]] &, 
                  alf]]]]];
   result = Merge[ 
    ParallelTable[
        Block[{tallies = Association[Map[# -> Association[] &, ops]],
               tt0 = SessionTime[], tt1, tmid = 0.0, tmidupdate = 0.0,
               tupdate = 0.0, tcm = 0.0}, 
       Do[Block[{x1 = latticeCoordinates[k1], x2, x22, xmiddle, x3,
                 w1, w2, w3, w1t, w2t, w3t, dx, face, f1, mid, flip},
        w1 = alf[[dir1, linearSiteIndex[x1]]];
        face = latticeDimensions;
        face[[dir1]] = 1;
        dy = vector[dir1, x1[[dir1]] - 1];
        (* Only use each pair of points once *)
        f1 = latticeIndex[x1 - dy, face];
        Do[
           x2 = latticeCoordinates[f2, face] + dy;
           w2 = alf[[dir1, linearSiteIndex[x2]]];
           dx = displace[x1 - x2];
           addTimeNull[tupdate, update2[tallies, "2t1", {Norm[dx]}, w1, w2]],
           {f2, f1}];
        Do[
           x2 = wrapIt[x1 + vector[dir1, i]];
           w2 = alf[[dir1, linearSiteIndex[x2]]];
           addTimeNull[tupdate, update2[tallies, "2l", {i}, w1, w2]];
           Do[
               If[dir1 != dir2,
                  Block[{dir3 = cross[dir1, dir2]},
                  w1t = alf[[dir2, linearSiteIndex[x1]]];
                  w2t = alf[[dir3, linearSiteIndex[x2]]];
                  {w1t, w2t} *= Signature[{dir1, dir2, dir3}];
                  addTimeNull[tupdate,
                              update2[tallies, "2t2", {i}, w1t, w2t]]]],
               {dir2, nd}];
           If[n==3,
              Do[
                  (*middle point*)
                  x3 = wrapIt[x1 + vector[dir1, j]];
                  w3 = alf[[dir1, linearSiteIndex[x3]]];
                  (*Reflection in dir1 direction*)
                  If[j < i - j, {w2, w3, w1} = -{w1, w3, w2}];
                  addTimeNull[tupdate,
                              update3[tallies, "3pl", {i, Max[j, i - j]},
                                      w1, w3, w2]],
                  {j, 0, i}];
              Do[If[dir2 != dir1, 
                    w1t = (alf[[dir2, linearSiteIndex[x1]]] + 
                          alf[[dir2, linearSiteIndex[shift[dir2, x1]]]])/2;
                    x22 = wrapIt[x1 + vector[dir1, If[flip, 1 - i, i]]];
                    w2t = alf[[dir2, linearSiteIndex[x22]]];
                    x3 = wrapIt[x1 + vector[dir1, If[flip, 1 - j, j]]];
                    w3t = alf[[dir2, linearSiteIndex[x3]]];
                    If[flip, w1t = -w1t];
                    addTimeNull[tupdate,
                                update3[tallies, "3pt1", {i - 1/2, j - 1/2},
                                        w1t, w3t, w2t]]],
                 {j, i}, {dir2, nd}, {flip, {False, True}}];
              Do[If[dir2 != dir1,
                    w1t = alf[[dir2, linearSiteIndex[x1]]];
                    w2t = alf[[dir2, linearSiteIndex[x2]]];
                    x3 = wrapIt[x1 + vector[dir1, j]];
                    w3 = (alf[[dir1, linearSiteIndex[x3]]] + 
                          alf[[dir1, linearSiteIndex[shift[dir2, x3]]]])/2;
                    (*Flip in dir1 direction*)
                    If[j + 1/2 < i - j - 1/2,
                       {w1t, w3, w2t} = {w2t, -w3, w1t}];
                    addTimeNull[
                        tupdate,
                        update3[tallies, "3pt2",
                                {i, Min[j + 1/2, i - j - 1/2]},
                                w1t, w3, w2t]]],
                 {j, 0, i - 1}, {dir2, nd}];
              face = latticeDimensions;
              face[[dir1]] = 1;
              dy = vector[dir1, x1[[dir1]] - 1 + Floor[i/2]];
              xmiddle = wrapIt[x1 + vector[dir1, Floor[i/2]]];
              Do[x3 = wrapIt[latticeCoordinates[f3, face] + dy];
                 w3 = If[Mod[i, 2] == 0, 
                         alf[[dir1, linearSiteIndex[x3]]],
                         (alf[[dir1, linearSiteIndex[x3]]] + 
                          alf[[dir1, linearSiteIndex[shift[dir1, x3]]]])/2];
                 dx = displace[xmiddle - x3];
                 addTimeNull[
                     tupdate,
                     update3[tallies, "3m1a", {i, Norm[dx]}, w1, w3, w2]],
                 {f3, latticeVolume[face]}];
              (* Ignore the "middle of the link" versus site. *)
              tt1 = SessionTime[];
              Do[If[dir1 != dir2,
                    w1t = alf[[d1, linearSiteIndex[x1]]];
                    w2t = alf[[d2, linearSiteIndex[x2]]];
                    x3 = wrapIt[x1 + vector[dir1, Floor[i/2]] +
                                vector[dir2, j]];
                    w3t = alf[[d3, linearSiteIndex[x3]]];
                    If[j < 0,
                       {w1t, w3t, w2t} *=
                       Map[If[# == dir2, -1, 1] &, {d1, d3, d2}]];
                    addTimeNull[
                        tcm,
                        {mid, flip} = componentMultiplet3[
                            dir1, dir2, d1, d3, d2]];
                    If[flip =!= Null,
                       {w2t, w3t, w1t} = {w1t, w3t, w2t}*flip];
                    addTimeNull[
                        tmidupdate,
                        update3[tallies, mid, {i, Abs[j]}, w1t, w3t, w2t]]],
                 {dir2, nd},
                 {j, -Floor[latticeDimensions[[dir2]]/2], 
                  latticeDimensions[[dir2]]/2},
                 {d1, nd}, {d2, nd}, {d3, nd}];
              tmid += SessionTime[] - tt1],
           {i, 0, latticeDimensions[[dir1]]/2}]],
          {k1, kernel, latticeVolume[], $KernelCount},
          {dir1, nd}];
       If[False, Print[Map[Map[Length[Flatten[#]]&, #]&, tallies]]];
       (* In principle, this is not needed, but
         the parallel kernel crashes for large lattices. *)
       If[distFlag,
          tallies = Map[Flatten, tallies, {2}]];
       If[OptionValue["timing"],
          Print["Timing:  ", {SessionTime[] - tt0, tmid - tcm - tmidupdate,
                              tupdate + tmidupdate, tcm}]];
       tallies],
      {kernel, $KernelCount}],
    Merge[#, If[distFlag, Flatten, Total]]&];
   Print["Total time: ", SessionTime[] - t0];
   Map[
       If[
           distFlag,
           Identity,
           {If[n == 2,
               (#[[2]] - #[[4]].#[[5]]/#[[1]])/(#[[1]] - 1),
               (#[[2]] + 2 ((fabc.#[[6]]).#[[5]]).#[[4]]/#[[1]]^2
                - Sum[#[[3 + i]].#[[6 + i]], {i, 3}]/#[[1]])/(#[[1]] - 1)],
            (* Use sample standard deviation (with Bessel's correction). *)
            #[[2]]/#[[1]], Sqrt[(#[[3]] - #[[2]]^2/#[[1]])/
                                (#[[1]]*(#[[1]]-1))]}&],
       result, {2}]]/;n==2||n==3;
