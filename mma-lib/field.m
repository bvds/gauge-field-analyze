opsList2 = {"2l", "2t1", "2t2", "2t11", "2t12", "2t13"};
opsList3 = {"3pl","3pt1", "3pt2", "3m1a", "3m1b", "3m2", 
            "3m3", "3m4", "3m5", "3m6", "3m7", "3m8", "3m9", "3m10"};
opsListAuto3 = {"3pl","3pt1", "3pt2", "3m1b", "3m2", 
            "3m4", "3m5", "3m6", "3m8", "3m9"};
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
 Block[{label1 = "\[CapitalDelta]",
        labelp = "\[CapitalDelta]", label2 = "r",
        label3 = "cov", labelm = "-cov"},
  Association[{
      "3pl" -> {1/2, -1, {label1, labelp, labelm}},
      "3pt1" -> {1, -1, {label1, labelp, labelm}},
      "3pt2" -> {1/2, 1, {label1, labelp, label3}},
      "3m1a" -> {Sqrt[1 + (#[[-2]]/#[[-1]])^2]&[latticeDimensions],
                 -1, {label1, label2, labelm}},
      "3m1b" -> {1, -1, {label1, label2, labelm}},
      "3m2" -> {1, -1, {label1, label2, labelm}}, 
      "3m3" -> {1, 1, {label1, label2, label3}},
      "3m4" -> {1, -1, {label1, label2, labelm}},
      "3m5" -> {1, -1, {label1, label2, labelm}},
      "3m6" -> {1, 1, {label1, label2, label3}},
      "3m7" -> {1, 1, {label1, label2, label3}},
      "3m8" -> {1, 1, {label1, label2, label3}},
      "3m9" -> {1, 1, {label1, label2, label3}},
      "3m10" -> {1, 1, {label1, label2, label3}}}]];


Options[linkCorrelators] = {"symmetric" -> False,
    "distFlag" -> True, "operators" -> All, "timing" -> False};
linkCorrelators[n_, OptionsPattern[]] :=
 Block[{
    distFlag = OptionValue["distFlag"],
    ops = Which[OptionValue["operators"] === All,
                If[n==2, opsList2, opsList3],
                OptionValue["operators"] === Automatic,
                If[n==2, opsList2, opsListAuto3],
                True,
                OptionValue["operators"]],
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
                 w1, w11, w2, w22, w3, dx, face, f1, mid, flip},
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
           addTimeNull[tupdate, update2[tallies, "2t1", {Norm[dx]}, w1, w2]];
           If[Length[opsList2] >= nd + 3,
              addTimeNull[tupdate, update2[tallies, opsList2[[3 + dir1]],
                                           {Norm[dx]}, w1, w2]]];
           If[False && Norm[dx]==Sqrt[13],
              Print["2t1 datum ", {dir1, x1, x2, dx}]],
           {f2, f1}];
        Do[
           x2 = wrapIt[x1 + vector[dir1, i]];
           w2 = alf[[dir1, linearSiteIndex[x2]]];
           addTimeNull[tupdate, update2[tallies, "2l", {i}, w1, w2]];
           Do[
               If[dir1 != dir2,
                  Block[{dir3 = cross[dir1, dir2]},
                  w11 = alf[[dir2, linearSiteIndex[x1]]];
                  w22 = alf[[dir3, linearSiteIndex[x2]]];
                  {w11, w22} *= Signature[{dir1, dir2, dir3}];
                  addTimeNull[tupdate,
                              update2[tallies, "2t2", {i}, w11, w22]]]],
               {dir2, nd}];
           If[n==3,
              Do[
                  (* Don't overwrite global value on reflection. *)
                  w11 = w1;
                  w22 = w2;
                  (*middle point*)
                  x3 = wrapIt[x1 + vector[dir1, j]];
                  w3 = alf[[dir1, linearSiteIndex[x3]]];
                  (*Reflection in dir1 direction*)
                  If[j > i - j, {w22, w3, w11} = -{w11, w3, w22}];
                  addTimeNull[tupdate,
                              update3[tallies, "3pl", {i, Min[j, i - j]},
                                      w11, w22, w3]],
                  {j, 0, i}];
              Do[If[dir2 != dir1, 
                    w11 = (alf[[dir1, linearSiteIndex[x1]]] + 
                          alf[[dir1, linearSiteIndex[shift[dir2, x1]]]])/2;
                    x22 = wrapIt[x1 + vector[dir1, If[flip, -i, i + 1]]];
                    w22 = alf[[dir2, linearSiteIndex[x22]]];
                    x3 = wrapIt[x1 + vector[dir1, If[flip, -j, j + 1]]];
                    w3 = alf[[dir2, linearSiteIndex[x3]]];
                    If[flip, w11 = -w11];
                    addTimeNull[tupdate,
                                update3[tallies, "3pt1", {i + 1/2, j + 1/2},
                                        w11, w22, w3]]],
                 {j, 0, i}, {dir2, nd}, {flip, {False, True}}];
              Do[If[dir2 != dir1,
                    w11 = alf[[dir2, linearSiteIndex[x1]]];
                    x22 = wrapIt[x1 + vector[dir1, i + 1]];
                    w22 = alf[[dir2, linearSiteIndex[x22]]];
                    x3 = wrapIt[x1 + vector[dir1, j]];
                    w3 = (alf[[dir1, linearSiteIndex[x3]]] + 
                          alf[[dir1, linearSiteIndex[shift[dir2, x3]]]])/2;
                    (*Flip in dir1 direction*)
                    If[j + 1/2 > i - j + 1/2,
                       {w11, w3, w22} = {w22, -w3, w11}];
                    addTimeNull[
                        tupdate,
                        update3[tallies, "3pt2",
                                {i + 1, Min[j + 1/2, i - j + 1/2]},
                                w11, w22, w3]]],
                 {j, 0, i}, {dir2, nd}];
              (* This is higher order in lattice size relative
                to other 3-point operators, so only go through
                the loop when needed. *)
              If[KeyExistsQ[tallies, "3m1a"],
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
                        update3[tallies, "3m1a", {i, Norm[dx]}, w1, w2, w3]],
                    {f3, latticeVolume[face]}]];
              (* Ignore the "middle of the link" versus site. *)
              tt1 = SessionTime[];
              Do[If[dir1 != dir2 &&
                    (* Only do j==0 once for longitudinal d3 *)
                    Not[j==0 && d3==dir1 && Mod[dir2 - dir1, 3] == 1],
                    w11 = alf[[d1, linearSiteIndex[x1]]];
                    w22 = alf[[d2, linearSiteIndex[x2]]];
                    x3 = wrapIt[x1 + vector[dir1, Floor[i/2]] +
                                vector[dir2, j]];
                    w3 = alf[[d3, linearSiteIndex[x3]]];
                    If[j < 0,
                       {w11, w3, w22} *=
                       Map[If[# == dir2, -1, 1] &, {d1, d3, d2}]];
                    addTimeNull[
                        tcm,
                        {mid, flip} = componentMultiplet3[
                            dir1, dir2, d1, d3, d2]];
                    If[flip =!= Null,
                       {w22, w3, w11} = {w11, w3, w22}*flip];
                    addTimeNull[
                        tmidupdate,
                        update3[tallies, mid, {i, Abs[j]}, w11, w22, w3]]],
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
   If[False,
      Print["3m1a for (2,0) ", result["3m1a"][{2,0}]];
      Print["3m1b for (2,0) ", result["3m1b"][{2,0}]]];
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

(* Use the fundamental representation *)
trace::usage = "Trace of the product of Hermitian matrices.  \
This allows for future optimization.";
trace[a_, b_] := Block[{x = Flatten[a].Flatten[Transpose[b]]},
                       If[Abs[Im[x]]>10^-8,
                          Print["non-hermetian arg", {a, b}];
                          Print[Column[Map[Short, Stack[_]]]];
                          Abort[]];
                       Re[x]];
trace[a__] := Tr[Dot[a]]/;Length[{a}]!=2;
opsList4 = {"t1", "m1", "t2", "m2", "t3"};
(* See arXiv:1207.0609v2 [hep-ph] 2 Oct 2012 for color multiplets. *)
colorList4 = {"1", "adjointS", "adjointA", "10S", "27", "0"};
linkCorrelators[4, OptionsPattern[]] :=
 Block[{
    distFlag = OptionValue["distFlag"],
    ops = Which[OptionValue["operators"] === All,
                Flatten[Outer[List, opsList4, colorList4], 1],
                OptionValue["operators"] === Automatic,
                Flatten[Outer[List, opsList4, colorList4], 1],
                True,
                OptionValue["operators"]],
    flf = Map[(-I SULog[#])&, gaugeField, {2}], 
    t0 = SessionTime[], t1, result,
    vector = Function[{dir, len}, Table[If[i == dir, len, 0], {i, nd}]],
    cross = Function[{k1, k2}, If[k1 != k2, 6 - k1 - k2, $Failed]],
    addHC = Function[x, x + ConjugateTranspose[x]],
    displace = Function[dx, 
       Mod[dx + latticeDimensions/2, latticeDimensions] -
       latticeDimensions/2],
    addTimeNull = If[
        OptionValue["timing"],
        Function[{timer, expr},
                 Block[{t1 = SessionTime[]},
                       expr; timer += SessionTime[] - t1],
                 {HoldAll, SequenceHold}],
        Null&],
    update},
   update = Function[
       {tallies, op, key, wa, wb, wc, wd},
       Block[{tab, tcd, tac, tbd, tad, tbc, sym, wab, wcd,
              p, d, name},
        tab = trace[wa, wb];
        tcd = trace[wc, wd];
        p["1"] = tab*tcd/(nc^2 -1);
        d["1"] = {wb*tcd, wa*tcd, wd*tab, wc*tab}/(nc^2 -1);
        Do[
            wab = wa.wb + sym*wb.wa;
            wcd = wc.wd + sym*wd.wc;
            name = If[sym > 0, "adjointS", "adjointA"];
            p[name] = trace[wab, wcd]*nc/(nc^2 - 2 - 2 sym);
            d[name] = {wb.wcd + sym*wcd.wb, wcd.wa + sym*wa.wcd,
                         wd.wab + sym*wab.wd, wab.wc + sym*wc.wab}*
                      nc/(nc^2 - 2 + 2 sym),
            {sym, {-1, 1}}];
        tac = trace[wa, wc];
        tbd = trace[wb, wd];
        tad = trace[wa, wd];
        tbc = trace[wb, wc];
        If[KeyExistsQ[tallies, {op, "10S"}],
           p["10S"] = (tac*tbd - tad*tbc)/2 - p["adjointA"];
           d["10S"] = {wc*tbd - wd*tbc, wd*tac - wc*tad,
                       wa*tbd - wb*tad, wb*tac - wa*tbc}/2 - d["adjointA"]];
        If[KeyExistsQ[tallies, {op, "27"}],
           p["27"] = (tac*tbd + tad*tbc +
                      2*Re[trace[wa, wc, wb, wd]])/4 -
                     p["adjointA"]*(nc+2)/(2 nc) - p["1"]*(nc+1)/(2 nc);
           d["27"] = {wc*tbd + wd*tbc + addHC[wc.wb.wd],
                      wd*tac + wc*tad + addHC[wd.wa.wc],
                      wa*tbd + wb*tad + addHC[wb.wd.wa],
                      wb*tac + wa*tbc + addHC[wa.wc.wb]}/4 -
                     d["adjointS"]*(nc+2)/(2 nc) - d["1"]*(nc+1)/(2 nc)];
        If[KeyExistsQ[tallies, {op, "0"}],
           p["0"] = (tac*tbd + tad*tbc -
                     2*Re[trace[wa, wc, wb, wd]])/4 -
                    p["adjointS"]*(nc+2)/(2 nc) - p["1"]*(nc+1)/(2 nc);
           d["0"] = {wc*tbd + wd*tbc - addHC[wc.wb.wd],
                     wd*tac + wc*tad - addHC[wd.wa.wc],
                     wa*tbd + wb*tad - addHC[wb.wd.wa],
                     wb*tac + wa*tbc - addHC[wa.wc.wb]}/4 -
                    d["adjointS"]*(nc+2)/(2 nc) - d["1"]*(nc+1)/(2 nc)];
        Do[
            If[KeyExistsQ[tallies, {op, colorm}],
               Block[{lt = Lookup[
                   tallies[{op, colorm}], Key[key],
                   If[distFlag,
                      {},
                      Join[{0, 0.0, 0.0},
                           Array[0.0&, {8, Length[wa], Length[wa]}]]]]},
                     If[Im[p[colorm]] != 0,
                        Print[{p[colorm], {op, key}}]];
                     Block[
                         {z = Join[{wa, wb, wc, wd}, d[colorm]]},
                         Do[
                             If[Norm[z[[i]] - ConjugateTranspose[z[[i]]]] > 10^-12,
                                Print["d not hermitian",
                                      {colorm, i, z[[i]]}]],
                             {i, Length[z]}]];
                     If[distFlag,
                        lt = {lt, p[colorm]},
                        lt += Join[{1, p[colorm], p[colorm]^2, wa, wb, wc, wd},
                                   d[colorm]]];
                     tallies[{op, colorm}][key] = lt]],
            {colorm, colorList4}]], {HoldFirst}];
   Print["ops: ", ops];
   If[False, Print[TableForm[
       Transpose[
           Map[Map[valueError[Mean[#], 
            StandardDeviation[#]/Sqrt[Length[#]]] &, Transpose[#]] &, 
                  flf]]]]];
   result = Merge[ 
    ParallelTable[
     Block[{tallies = Association[Map[# -> Association[] &, ops]],
               tt0 = SessionTime[], tt1, tupdate = 0.0}, 
      Do[Block[{x1 = latticeCoordinates[k1], x2, dx, face, f1, w1, dy, w2},
        w1 = flf[[dir1, linearSiteIndex[x1]]];
        face = latticeDimensions;
        face[[dir1]] = 1;
        dy = vector[dir1, x1[[dir1]] - 1];
        (* Only use each pair of points once *)
        f1 = latticeIndex[x1 - dy, face];
        Do[
            x2 = latticeCoordinates[f2, face] + dy;
            w2 = flf[[dir1, linearSiteIndex[x2]]];
            dx = displace[x1 - x2];
            addTimeNull[tupdate, update[tallies, "t2", {Norm[dx]},
                                        w1, w1, w2, w2]],
            {f2, f1}];
        Do[
            x2 = wrapIt[x1 + vector[dir1, i]];
            w2 = flf[[dir1, linearSiteIndex[x2]]];
            addTimeNull[tupdate, update[tallies, "l", {i}, w1, w1, w2, w2]];
            Do[
                If[dir1 != dir2,
                 Block[
                   {dir3 = cross[dir1, dir2], w12, w22, w13, w23, sig},
                  sig = Signature[{dir1, dir2, dir3}];
                  w12 = flf[[dir2, linearSiteIndex[x1]]];
                  w22 = flf[[dir2, linearSiteIndex[x2]]];
                  w13 = flf[[dir3, linearSiteIndex[x1]]];
                  w23 = flf[[dir3, linearSiteIndex[x2]]];
                  addTimeNull[tupdate,
                              update[tallies, "t1", {i},
                                     w12, sig*w13, w22, sig*w23]];
                  addTimeNull[tupdate,
                              update[tallies, "m1", {i},
                                     w12, w1, w22, w2]];
                  addTimeNull[tupdate,
                              update[tallies, "m2", {i},
                                     w1, w1, w22, w22]];
                  addTimeNull[tupdate,
                              update[tallies, "m2", {i},
                                     -w12, -w2, w12, w12]];
                  addTimeNull[tupdate,
                              update[tallies, "t3", {i},
                                     sig*w13, sig*w13, w22, w22]]]],
                {dir2, nd}],
            {i, 0, latticeDimensions[[dir1]]/2}]],
         {k1, kernel, latticeVolume[], $KernelCount},
         {dir1, nd}];
       If[False, Print[Map[Map[Length[Flatten[#]]&, #]&, tallies]]];
       (* In principle, this is not needed, but
         the parallel kernel crashes for large lattices. *)
       If[distFlag,
          tallies = Map[Flatten, tallies, {2}]];
       If[OptionValue["timing"],
          Print["Timing:  ", {SessionTime[] - tt0, tupdate}]];
       tallies],
     {kernel, $KernelCount}],
    Merge[#, If[distFlag, Flatten, Total]]&];
   Print["Total time: ", SessionTime[] - t0];
   If[False, Print["Result:  ", result]];
   Map[
       If[
           distFlag,
           Identity,
           (* Drop terms with more than one average value. *)
           {(#[[2]] -
             Sum[trace[#[[3 + i]], #[[7 + i]]],
                 {i, 4}]/#[[1]])/(#[[1]] - 1),
            (* Use sample standard deviation (with Bessel's correction). *)
            If[Im[#[[2]]] != 0, Print["bad ", #]];
            #[[2]]/#[[1]], Sqrt[(#[[3]] - #[[2]]^2/#[[1]])/
                                (#[[1]]*(#[[1]]-1))]}&],
       result, {2}]];
