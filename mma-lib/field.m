opsList2 = {"2l", "2t1", "2t2", "2t11", "2t12", "2t13"};
opsList3 = {"3pl","3pt1", "3pt2", "3m1a", "3m1b", "3m2", 
            "3m3", "3m4", "3m5", "3m6", "3m7", "3m8", "3m9", "3m10"};
componentMultiplet3 = {
         {"3m1b", {"l", "l", "l"}, Null},
         {"3m2", {"l", "t", "l"}, Null},
         {"3m3", {"l", "l", "t"}, Null},
         {"3m4", {"l", "t", "t"}, Null},
         {"3m5", {"l", "z", "z"}, Null},
         {"3m2", {"t", "l", "l"}, {1, -1, -1}},
         {"3m6", {"t", "t", "l"}, Null},
         {"3m4", {"t", "l", "t"}, {1, -1, 1}},
         {"3m7", {"t", "t", "t"}, Null},
         {"3m8", {"t", "z", "z"}, Null},
         {"3m9", {"z", "z", "l"}, Null},
         {"3m10", {"z", "z", "t"}, Null},
         {"3m5", {"z", "l", "z"}, {1, -1, 1}},
         {"3m8", {"z", "t", "z"}, {1, 1, 1}}};

linkCorrPlot3::usage = "Helper for plot routines, returning:  aspect ratio, sign of the data, and axis labels.";
linkCorrPlot3[] :=
 Block[{label1 = "\[CapitalDelta]",
        labelp = "\[Delta]", label2 = "r",
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


Options[linkCorrelators] = {
    (* See Malin Sjödahl, arXiv:0906.1121v2 [hep-ph] 31 Jul 2009
      for a discussion why these are zero.  I believe that this
      is because the vacuum is even under charge conjugation.*)
    "anomoly" -> False,
    "distFlag" -> False, "timing" -> False,
    (* Drop any that are zero under symmetries *)
    "operators" -> Automatic};
linkCorrelators[n_, OptionsPattern[]] :=
 Block[{
    distFlag = OptionValue["distFlag"],
    ops = Which[OptionValue["operators"] === All,
                If[n==2, opsList2, opsList3],
                OptionValue["operators"] === Automatic,
                If[n==2,
                   opsList2,
                   (* Remove cases that are zero. *)
                   Complement[opsList3,
                              If[OptionValue["anomoly"],
                                 {"3m1a", "3m1b", "3m6", "3m9"},
                                 (* remove "3m1a" because it is big *)
                                 {"3m1a", "3m3", "3m7", "3m10"}]]],
                True,
                OptionValue["operators"]],
    (* Subtract average value separately in each direction *)
    alf = Map[Block[{x=Mean[#]}, Map[(#-x)&, #]]&,
              Map[2*Im[Map[Tr, SUGenerators[].SULog[#]]] &, gaugeField, {2}]], 
    getf = Function[{flf, dir, x},
                      flf[[dir, linearSiteIndex[x]]], {HoldFirst}],
    t0 = SessionTime[], t1, result,
    half = IdentityMatrix[nd]/2, 
    vector = Function[{dir, len}, Table[If[i == dir, len, 0], {i, nd}]],
    cross = Function[{k1, k2}, If[k1 != k2, 6 - k1 - k2, $Failed]],
    displace = Function[dx, 
       Mod[dx + latticeDimensions/2, latticeDimensions] -
       latticeDimensions/2],
    fabc = N[If[OptionValue["anomoly"], SUSymmetric[], SUStructure[]]],
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
                                {0, 0.0, 0.0}]],
                w12 = w1.w2},
                If[distFlag,
                   lt = {lt, w12},
                   lt += {1, w12, w12^2}];
                tallies[op][key] = lt]], {HoldFirst}];
   update3 = Function[
       {tallies, op, key, w1, w2, w3}, 
       If[KeyExistsQ[tallies, op], 
          (* Right multiplication is faster for SparseArray. *)
          Block[{f321 = ((fabc.w3).w2).w1, lt = Lookup[
              tallies[op], Key[key],
              If[distFlag, {},
                 {0, 0.0, 0.0}]]},
                If[distFlag,
                   lt = {lt, f321},
                   lt += {1, f321, f321^2}];
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
       Do[Block[{x1 = latticeCoordinates[k1], dx, dy, face, f1, mid, flip,
                 x0, x2, x22, x3,
                 w1 = Association[], w2 =Association[], w3 = Association[],
                 wa, wb, wc},
        w1["l"] = getf[alf, dir1, x1];
        face = latticeDimensions;
        face[[dir1]] = 1;
        dy = vector[dir1, x1[[dir1]] - 1];
        (* Only use each pair of points once *)
        f1 = latticeIndex[x1 - dy, face];
        Do[
           x2 = latticeCoordinates[f2, face] + dy;
           w2["l"] = getf[alf, dir1, x2];
           dx = displace[x1 - x2];
           addTimeNull[tupdate,
                       update2[tallies, "2t1", {Norm[dx]}, w1["l"], w2["l"]]];
           If[Length[opsList2] >= nd + 3,
              addTimeNull[tupdate, update2[tallies, opsList2[[3 + dir1]],
                                           {Norm[dx]}, w1["l"], w2["l"]]]];
           If[False && Norm[dx]==Sqrt[13],
              Print["2t1 datum ", {dir1, x1, x2, dx}]],
           {f2, f1}];
        Do[
           x2 = shift[dir1, x1, i];
           w2["l"] = getf[alf, dir1, x2];
           addTimeNull[tupdate,
                       update2[tallies, "2l", {i}, w1["l"], w2["l"]]];
           Do[
               If[dir1 != dir2,
                 Block[{dir3 = cross[dir1, dir2]},
                  wa = getf[alf, dir2, x1];
                  wb = getf[alf, dir3, x2];
                  {wa, wb} *= Signature[{dir1, dir2, dir3}];
                  addTimeNull[tupdate,
                              update2[tallies, "2t2", {i}, wa, wb]]]],
               {dir2, nd}];
           If[n==3,
              Do[
                  (* Don't overwrite global value on reflection. *)
                  wa = w1["l"];
                  wb = w2["l"];
                  (* middle point *)
                  x3 = shift[dir1, x1, j];
                  wc = getf[alf, dir1, x3];
                  (* Reflection in dir1 direction *)
                  If[j > i - j, {wb, wa, wc} = -{wa, wb, wc}];
                  addTimeNull[tupdate,
                              update3[tallies, "3pl", {i, Min[j, i - j]},
                                      wa, wb, wc]],
                  {j, 0, i}];
              Do[If[dir2 != dir1, 
                    wa = (getf[alf, dir1, x1] +
                          getf[alf, dir1, shift[dir2, x1]])/2;
                    x22 = shift[dir1, x1, If[flip, -i, i + 1]];
                    wb = getf[alf, dir2, x22];
                    x3 = shift[dir1, x1, If[flip, -j, j + 1]];
                    wc = getf[alf, dir2, x3];
                    If[flip, wa = -wa];
                    addTimeNull[tupdate,
                                update3[tallies, "3pt1", {i + 1/2, j + 1/2},
                                        wa, wb, wc]]],
                 {j, 0, i}, {dir2, nd}, {flip, {False, True}}];
              Do[If[dir2 != dir1,
                    wa = getf[alf, dir2, x1];
                    x22 = shift[dir1, x1, i + 1];
                    wb = getf[alf, dir2, x22];
                    x3 = shift[dir1, x1, j];
                    wc = (getf[alf, dir1, x3] + 
                          getf[alf, dir1, shift[dir2, x3]])/2;
                    (* Flip in dir1 direction *)
                    If[j + 1/2 > i - j + 1/2,
                       {wb, wa, wc} = {wa, wb, -wc}];
                    addTimeNull[
                        tupdate,
                        update3[tallies, "3pt2",
                                {i + 1, Min[j + 1/2, i - j + 1/2]},
                                wa, wb, wc]]],
                 {j, 0, i}, {dir2, nd}];
              (* This is higher order in lattice size relative
                to other 3-point operators, so only go through
                the loop when needed.  *)
              If[KeyExistsQ[tallies, "3m1a"],
                 face = latticeDimensions;
                 face[[dir1]] = 1;
                 dy = vector[dir1, x1[[dir1]] - 1 + Floor[i/2]];
                 x0 = shift[dir1, x1, Floor[i/2]];
                 Do[x3 = wrapIt[latticeCoordinates[f3, face] + dy];
                    wc = If[Mod[i, 2] == 0, 
                            getf[alf, dir1, x3],
                            (getf[alf, dir1, x3] + 
                             getf[alf, dir1, shift[dir1, x3]])/2];
                    dx = displace[x0 - x3];
                    addTimeNull[
                        tupdate,
                        update3[tallies, "3m1a", {i, Norm[dx]},
                                w1["l"], w2["l"], wc]],
                    {f3, latticeVolume[face]}]];
              (* Ignore the "middle of the link" versus site. *)
              tt1 = SessionTime[];
              Do[If[dir1 != dir2,
                    Block[{dir3 = cross[dir1, dir2]},
                     w1["t"] = getf[alf, dir2, x1];
                     w1["z"] = getf[alf, dir3, x1];
                     w2["t"] = getf[alf, dir2, x2];
                     w2["z"] = getf[alf, dir3, x2];
                     x3 = wrapIt[x1 + vector[dir1, Floor[i/2]] +
                                 vector[dir2, j]];
                     w3["l"] = getf[alf, dir1, x3];
                     w3["t"] = getf[alf, dir2, x3];
                     w3["z"] = getf[alf, dir3, x3];
                     If[j < 0,
                        w1["t"] *= -1;
                        w2["t"] *= -1;
                        w3["t"] *= -1];
                     addTimeNull[
                         tcm,
                         Scan[Apply[Function[
                             {mid, dirs, flip},
                             If[
                                 KeyExistsQ[tallies, mid] &&
                                 (* Only do j==0 once for longitudinal w3 *)
                                 j!=0 || dirs[[3]] != "l" ||
                                 Signature[{dir1, dir2, dir3}] > 0,
                                 {wa, wb, wc} = MapThread[
                                     #1[#2]&, {{w1, w2, w3}, dirs}];
                                 If[flip =!= Null,
                                    {wb, wa, wc} = {wa, wb, wc}*flip];
                                 addTimeNull[
                                     tmidupdate,
                                     update3[tallies, mid, {i, Abs[j]},
                                             wa, wb, wc]]]],
                                    #]&,
                                   componentMultiplet3]
                     ]]],
                 {dir2, nd},
                 {j, -Floor[latticeDimensions[[dir2]]/2], 
                  latticeDimensions[[dir2]]/2}];
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
            (* Use sample standard deviation (with Bessel's correction). *)
            {#[[2]]/#[[1]], Sqrt[(#[[3]] - #[[2]]^2/#[[1]])/
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
opsList22 = {"t1", "m1", "l", "t2", "m2", "t3"};
(* See arXiv:1207.0609v2 [hep-ph] 2 Oct 2012 for color multiplets.
   arXiv:hep-ph/9803241v1 4 Mar 1998 gives more complete discussion. *)
colorList4[] := colorList4[nc];
colorList4[nc_] := {"1", "8S", "8A", "Re10", "27",
                    If[nc > 3, "0", Nothing]};
colorList4A = {"A1", "A2", "A3"};
SetAttributes[updateColor4, HoldFirst];
updateColor4[tallies_, op_, key_, distFlag_, wa_, wb_, wc_, wd_] :=
 Block[{tab, tcd, tac, tbd, tad, tbc, p = Association[]},
  tab = trace[wa, wb];
  tcd = trace[wc, wd];
  p["1"] = 4*tab*tcd/(nc^2 -1);
  p["8S"] = 2*(2*Re[trace[wa, wb, wc, wd]] +
                     2*Re[trace[wb, wa, wc, wd]] - 4*tab*tcd/nc)*
                  nc/(nc^2 - 4);
  p["8A"] = 2*(2*Re[trace[wa, wb, wc, wd]] -
                     2*Re[trace[wb, wa, wc, wd]])/nc;
  tac = trace[wa, wc];
  tbd = trace[wb, wd];
  tad = trace[wa, wd];
  tbc = trace[wb, wc];
  If[KeyExistsQ[tallies, {op, "Re10"}],
     (* Taking the real part of the 10 *)
     p["Re10"] = tac*tbd - tad*tbc - p["8A"]/2];
  If[KeyExistsQ[tallies, {op, "27"}],
     p["27"] = tac*tbd + tad*tbc +
               2*Re[trace[wa, wc, wb, wd]] -
               p["8S"]*(nc-2)/(2 nc) - p["1"]*(nc-1)/(2 nc)];
  If[KeyExistsQ[tallies, {op, "0"}],
     p["0"] = tac*tbd + tad*tbc -
              2*Re[trace[wa, wc, wb, wd]] -
              p["8S"]*(nc+2)/(2 nc) - p["1"]*(nc+1)/(2 nc)];
  If[KeyExistsQ[tallies, {op, "A1"}],
     p["A1"] = 2*Im[trace[wa, wb, wc, wd]]];
  If[KeyExistsQ[tallies, {op, "A2"}],
     p["A2"] = 2*Im[trace[wa, wb, wd, wc]]];
  If[KeyExistsQ[tallies, {op, "A3"}],
     p["A3"] = 2*Im[trace[wa, wc, wb, wd]]];
  Do[
      If[KeyExistsQ[tallies, {op, colorm}],
         Block[{lt = Lookup[
             tallies[{op, colorm}], Key[key],
             If[distFlag, {}, {0, 0.0, 0.0}]]},
               If[Im[p[colorm]] != 0,
                  Print[{p[colorm], {op, key}}]];
               If[distFlag,
                  lt = {lt, p[colorm]},
                  lt += {1, p[colorm], p[colorm]^2}];
               tallies[{op, colorm}][key] = lt]],
      (* List of color multiplets defined for p. *)
      {colorm, Keys[p]}]];

linkCorrelators[22, OptionsPattern[]] :=
 Block[{
    distFlag = OptionValue["distFlag"],
    ops = Which[OptionValue["operators"] === All,
                Flatten[Outer[List, opsList22,
                              If[OptionValue["anomoly"],
                                 colorList4A, colorList4[]]], 1],
                OptionValue["operators"] === Automatic,
                (* Remove miltiplets that are zero. *)
                Complement[
                    Flatten[Outer[List, opsList22,
                                  If[OptionValue["anomoly"],
                                     colorList4A, colorList4[]]], 1],
                    Flatten[Outer[List,
                                  {"l", "t2", "m2", "t3"},
                                  {"8A", "Re10"}], 1]],
                True,
                OptionValue["operators"]],
    (* Subtract average value separately in each direction *)
    flf = Map[Block[{x=Mean[#]}, Map[(#-x)&, #]]&,
                   Map[(-I SULog[#])&, gaugeField, {2}]],
    getf = Function[{flf, dir, x},
                      flf[[dir, linearSiteIndex[x]]], {HoldFirst}],
    t0 = SessionTime[], t1, result,
    vector = Function[{dir, len}, Table[If[i == dir, len, 0], {i, nd}]],
    cross = Function[{k1, k2}, If[k1 != k2, 6 - k1 - k2, $Failed]],
    displace = Function[dx, 
       Mod[dx + latticeDimensions/2, latticeDimensions] -
       latticeDimensions/2],
    addTimeNull = If[
        OptionValue["timing"],
        Function[{timer, expr},
                 Block[{t1 = SessionTime[]},
                       expr; timer += SessionTime[] - t1],
                 {HoldAll, SequenceHold}],
        Null&]},
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
      Do[Block[{x1 = latticeCoordinates[k1], x2, dx, face, f1, dy,
                     w1l, w2l},
        w1l = getf[flf, dir1, x1];
        face = latticeDimensions;
        face[[dir1]] = 1;
        dy = vector[dir1, x1[[dir1]] - 1];
        (* Only use each pair of points once *)
        f1 = latticeIndex[x1 - dy, face];
        Do[
            x2 = latticeCoordinates[f2, face] + dy;
            w2l = getf[flf, dir1, x2];
            dx = displace[x1 - x2];
            addTimeNull[tupdate, updateColor4[tallies, "t2", {Norm[dx]},
                                        distFlag, w1l, w1l, w2l, w2l]],
            {f2, f1}];
        Do[
            x2 = shift[dir1, x1, i];
            w2l = getf[flf, dir1, x2];
            addTimeNull[tupdate, updateColor4[tallies, "l", {i}, distFlag,
                                              w1l, w1l, w2l, w2l]];
            Do[
                If[dir1 != dir2,
                 Block[
                     {dir3 = cross[dir1, dir2], w1t, w2t, w1z, w2z,
                      w1tt, w2tt, w1zz, w2zz, sig},
                  sig = Signature[{dir1, dir2, dir3}];
                  w1t = getf[flf, dir2, x1];
                  w2t = getf[flf, dir2, x2];
                  w1tt = getf[flf, dir2, shift[dir2, x1, -1]];
                  w2tt = getf[flf, dir2, shift[dir2, x2, -1]];
                  If[sig > 0,
                     w1z = getf[flf, dir3, x1];
                     w2z = getf[flf, dir3, x2];
                     w1zz = getf[flf, dir3, shift[dir3, x1, -1]];
                     w2zz = getf[flf, dir3, shift[dir3, x2, -1]];
                     addTimeNull[tupdate,
                                 updateColor4[tallies, "t1", {i}, distFlag,
                                              w1t, w1z, w2t, w2z]];
                     addTimeNull[tupdate,
                                 updateColor4[tallies, "t1", {i}, distFlag,
                                              w1z, -w1tt, w2z, -w2tt]];
                     addTimeNull[tupdate,
                                 updateColor4[tallies, "t1", {i}, distFlag,
                                              -w1zz, w1t, -w2zz, w2t]];
                     addTimeNull[tupdate,
                                 updateColor4[tallies, "t1", {i}, distFlag,
                                              -w1tt, -w1zz, -w2tt, -w2zz]];
                     addTimeNull[tupdate,
                                 updateColor4[tallies, "t3", {i}, distFlag,
                                              w1z, w1zz, w2t, w2tt]];
                     addTimeNull[tupdate,
                                 updateColor4[tallies, "t3", {i}, distFlag,
                                              w2z, w2zz, w1t, w1tt]]];
                  addTimeNull[tupdate,
                              updateColor4[tallies, "m1", {i}, distFlag,
                                     w1t, w1l, w2t, w2l]];
                  addTimeNull[tupdate,
                              updateColor4[tallies, "m1", {i}, distFlag,
                                     -w1tt, w1l, -w2tt, w2l]];
                  If[i > 0,
                     addTimeNull[tupdate,
                                 updateColor4[tallies, "m2", {i-1/2}, distFlag,
                                              w1l, w1l, w2t, w2tt]]];
                  addTimeNull[tupdate,
                              updateColor4[tallies, "m2", {i+1/2}, distFlag,
                                           -w2l, -w2l, w1t, w1tt]]]],
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
            (* Use sample standard deviation (with Bessel's correction). *)
            {#[[2]]/#[[1]], Sqrt[(#[[3]] - #[[2]]^2/#[[1]])/
                                (#[[1]]*(#[[1]]-1))]}&],
       result, {2}]];

opsList4x = {"4hr", "4hh", "4hz", "4zr", "4zz"};
linkCorrelators["4x", OptionsPattern[]] :=
 Block[{
    distFlag = OptionValue["distFlag"],
    ops = Which[OptionValue["operators"] === All,
                Flatten[Outer[List, opsList4x,
                              If[OptionValue["anomoly"],
                                 colorList4A, colorList4[]]], 1],
                OptionValue["operators"] === Automatic,
                Flatten[Outer[List, opsList4x,
                              If[OptionValue["anomoly"],
                                 colorList4A,
                                 (* Must be symmetric under w1 <-> w2 *)
                                 Complement[colorList4[],
                                            {"8A", "Re10"}]]], 1],
                True,
                OptionValue["operators"]],
    (* Subtract average value separately in each direction *)
    flf = Map[Block[{x=Mean[#]}, Map[(#-x)&, #]]&,
              Map[(-I SULog[#])&, gaugeField, {2}]], 
    getf = Function[{flf, dir, x},
                      flf[[dir, linearSiteIndex[x]]], {HoldFirst}],
    t0 = SessionTime[], t1, result,
    cross = Function[{k1, k2}, If[k1 != k2, 6 - k1 - k2, $Failed]],
    displace = Function[dx, 
       Mod[dx + latticeDimensions/2, latticeDimensions] -
       latticeDimensions/2],
    addTimeNull = If[
        OptionValue["timing"],
        Function[{timer, expr},
                 Block[{t1 = SessionTime[]},
                       expr; timer += SessionTime[] - t1],
                 {HoldAll, SequenceHold}],
        Null&]},
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
      Do[
       Block[{x0 = latticeCoordinates[k1]},
        Do[If[dir1 != dir2,
         Block[{dir3 = cross[dir1, dir2], sig,
                x1, x11, x2, x3, x33, x4, x44,
                w1h, w1z, w1zz, w11h, w2h, w2z, w2zz,
                w3r, w3z, w33h, w33z, w33zz, w4h, w4r, w4z, w4zz, w44r},
          sig = Signature[{dir1, dir2, dir3}];
          x1 = shift[dir1, x0, i];
          w1h = getf[flf, dir2, x1];
          w1z = getf[flf, dir3, x1];
          w1zz = (w1z + getf[flf, dir3, shift[dir3, x1, -1]])/2;
          x11 = shift[dir1, x0, i+1];
          w11h = getf[flf, dir2, x11];
          x2 = shift[dir1, x0, -i];
          w2h = getf[flf, dir2, x2];
          w2z = getf[flf, dir3, x2];
          w2zz = (w2z + getf[flf, dir3, shift[dir3, x2, -1]])/2;
          Do[
              x3 = shift[dir2, x0, j];
              w3r = getf[flf, dir2, x3];
              w3z = getf[flf, dir3, x3];
              x33 = shift[dir2, x0, j+1];
              w33h = getf[flf, dir1, x33];
              w33z = getf[flf, dir3, x33];
              w33zz = (w3z + getf[flf, dir3, shift[dir3, x33, -1]])/2;
              x4 = shift[dir2, x0, -j];
              w4r = getf[flf, dir2, x4];
              w4h = getf[flf, dir1, x4];
              w4z = getf[flf, dir3, x4];
              w4zz = (w4z + getf[flf, dir3, shift[dir3, x4, -1]])/2;
              x44 = shift[dir2, x0, -j-1];
              w44r = getf[flf, dir2, x44];
              addTimeNull[tupdate,
                          updateColor4[tallies, "4hr", {2*i, 2*j}, distFlag,
                                       w1h, w2h, w3r, w4r]];
              addTimeNull[tupdate,
                          updateColor4[tallies, "4hz", {2*i, 2*j+1}, distFlag,
                                       w1h, w2h, sig*w33zz, sig*w4zz]];
              addTimeNull[tupdate,
                          updateColor4[tallies, "4zr", {2*i, 2*j+1}, distFlag,
                                       sig*w1zz, sig*w2zz, w3r, w44r]];
              If[sig>0,
                 If[j <= i,
                 addTimeNull[tupdate,
                             updateColor4[tallies, "4hh", {2*i+1, 2*j+1},
                                          distFlag,
                                          w11h, w2h, w33h, w4h]],
                    addTimeNull[tupdate,
                                updateColor4[tallies, "4hh", {2*j+1, 2*i+1},
                                             distFlag,
                                             w33h, w4h, w11h, w2h]]];
                 If[j <= i,
                    addTimeNull[tupdate,
                                updateColor4[tallies, "4zz", {2*i, 2*j},
                                             distFlag,
                                             w1z, w2z, w3z, w4z]],
                    addTimeNull[tupdate,
                                updateColor4[tallies, "4zz", {2*j, 2*i},
                                             distFlag,
                                             w3z, w4z, w1z, w2z]]]],
              {j, 0, latticeDimensions[[dir1]]/4}]]],
           {dir1, nd}, {i, 0, latticeDimensions[[dir1]]/4},
           {dir2, nd}]],
       {k1, kernel, latticeVolume[], $KernelCount}];
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
            (* Use sample standard deviation (with Bessel's correction). *)
            {#[[2]]/#[[1]], Sqrt[(#[[3]] - #[[2]]^2/#[[1]])/
                                (#[[1]]*(#[[1]]-1))]}&],
       result, {2}]];
