makeFieldCorrelators[] :=
 Block[{agf = 
        Map[2*Im[Map[Tr, SUGenerators[].SULog[#]]] &, gaugeField, {2}], 
        t0 = SessionTime[], t1, result},
  Print[Dimensions[agf]]; 
  Print[TableForm[
      Transpose[
          Map[Map[valueError[Mean[#], 
          StandardDeviation[#]/Sqrt[Length[#]]] &, Transpose[#]] &, 
      agf]]]]; 
  result = Map[Merge[#, Total] &, 
    Transpose[
     ParallelTable[
      Block[{tallies1 = Association[], tallies2 = Association[], 
        half = IdentityMatrix[nd]/2}, 
       Do[Block[{x1 = latticeCoordinates[k1], 
                 x2 = latticeCoordinates[k2], w1, w2, dx, lt, key},
         (* measure displacement relative to the middle of the link *)
         dx = Mod[(x1 + half[[dir1]]) - (x2 + half[[dir2]]) + 
                  latticeDimensions/2, latticeDimensions] - 
              latticeDimensions/2;
         w1 = agf[[dir1, linearSiteIndex[x1]]];
         w2 = agf[[dir2, linearSiteIndex[x2]]]; 
         If[dx[[dir1]] < 0, w1 = -w1];
         If[dx[[dir2]] < 0, w2 = -w2]; 
         key = If[
           dir1 == dir2, {Abs[dx[[dir1]]], 
            Sqrt[Total[dx^2] - dx[[dir1]]^2]}, 
           Block[{dxplane = Abs[dx[[{dir1, dir2}]]]}, {Max[dxplane], 
             Min[dxplane], Sqrt[Total[dx^2] - Total[dxplane^2]]}]]; 
         lt = Lookup[
             If[dir1 == dir2, tallies1, tallies2], 
             Key[key],
             {0, Table[0.0, {nc^2 - 1}, {nc^2 - 1}],
              Table[0.0, {nc^2 - 1}], Table[0.0, {nc^2 - 1}]}]; 
         lt += {1, Outer[Times, w1, w2], w1, w2}; 
         If[dir1 == dir2, tallies1[key] = lt, 
            tallies2[key] = lt]],
          {k1, kernel, latticeVolume[], $KernelCount},
          {k2, latticeVolume[]}, {dir1, nd}, {dir2, nd}];
       {tallies1, tallies2}], {kernel, $KernelCount}]]];
  Print["Time: ", SessionTime[] - t0];
  Map[(#[[2]] - Outer[Times, #[[3]], #[[4]]]/#[[1]])/(#[[1]] - 1) &, 
  result, {2}]];
