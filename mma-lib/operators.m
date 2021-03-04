plaquette[dir1_, dir2_, coords_List] := Tr[
    getLink[dir1, coords].
    getLink[dir2, shift[dir1, coords]].
    ConjugateTranspose[getLink[dir1, shift[dir2, coords]]].
    ConjugateTranspose[getLink[dir2, coords]]];
averagePlaquette[] :=
    Sum[Re[plaquette[dir1, dir2, latticeCoordinates[k]]],
	     {k, latticeVolume[]}, {dir1, 2, nd},
	     {dir2, dir1 - 1}]*2/(latticeVolume[]*nc*nd*(nd - 1));
averagePlaquette[dir1_, dir2_] :=
    Sum[Re[plaquette[dir1, dir2, latticeCoordinates[k]]],
	     {k, latticeVolume[]}]/(latticeVolume[]*nc);


smearLinks::usage = "For a given lattice, smear spatial link fields using Eqn. (23) in B. Lucini, M. Teper, U. Wenger https://arxiv.org/abs/hep-lat/0404008.  Returns the smeared lattice configuration.";
smearLinks::elbow = "Elbow operators not implemented.";
timeSmearLinks = 0;
smearLinks[dir0_, pa_:0.30, pd_:0.12] :=
    Block[
        {gf = gaugeField, i=1, t0 = SessionTime[], result},
        Which[
            True,
            (* Non-parallel version with sum over staples fixed *)
            Do[
                If[dir1 != dir0, Block[
                    {coord = latticeCoordinates[k], ff, xx, t0, z},
                    ff = Sum[
                        xx = shift[dir2, coord, -1];
                        getLink[dir1, coord] + pa*(
                            getLink[dir2, coord] .
                            getLink[dir1, shift[dir2, coord]] .
                            ConjugateTranspose[getLink[dir2, shift[dir1, coord]]] +
                            ConjugateTranspose[getLink[dir2, xx]] .
                            getLink[dir1, xx].
                            getLink[dir2, shift[dir1, xx]]),
                        {dir2, Complement[Range[nd], {dir0, dir1}]}];
                    If[nd > 3,
                       Message[smearLinks::elbow];
                       Return[$Failed]];
                    If[coord == {1,1,1} && i++ <5 && False,
                       Print[{{dir0, dir1, dir2}, ff}]];
                    t0 = SessionTime[];
                    gf[[dir1, linearSiteIndex[coord]]] =
                    SUStapleMinimum[
                        -ConjugateTranspose[ff],
                        printLevel -> If[coord == {1, 1, 1} && False, 1, 0]];
                    tMinimum += SessionTime[] - t0;
                    {dir1, coord, z}]],
                {k, latticeVolume[]}, {dir1, nd}];
            timeSmearLinks += SessionTime[] - t0,
            True,
            (* Parallel version *)
            result = ParallelTable[
                If[dir1 != dir0, Block[
                    {coord = latticeCoordinates[k], ff, xx, t0, z},
                    ff = Sum[
                        xx = shift[dir2, coord, -1];
                        getLink[dir1, coord] + pa*(
                            getLink[dir2, coord] .
                            getLink[dir1, shift[dir2, coord]] .
                            ConjugateTranspose[getLink[dir2, shift[dir1, coord]]] +
                            ConjugateTranspose[getLink[dir2, xx]] .
                            getLink[dir1, xx].
                            getLink[dir2, shift[dir1, xx]]),
                        {dir2, Complement[Range[nd], {dir0, dir1}]}];
                    If[nd > 3,
                       Message[smearLinks::elbow];
                       Return[$Failed]];
                    If[coord == {1,1,1} && i++ <5 && False,
                       Print[{{dir0, dir1, dir2}, ff}]];
                    t0 = SessionTime[];
                    z = SUStapleMinimum[
                        -ConjugateTranspose[ff],
                        printLevel -> If[coord == {1, 1, 1} && False, 1, 0]];
                    tMinimum += SessionTime[] - t0;
                    {dir1, linearSiteIndex[coord], z}],
                   Nothing],
                {k, latticeVolume[]}, {dir1, nd}];
            Scan[(gf[[#[[1]], #[[2]]]] = #[[3]])&,
                                 result, {2}];
            timeSmearLinks += SessionTime[] - t0;
        ];
        gf];

simpleBlock::usage = "For a given lattice, create a new lattice that is blocked in the spatial directions, using Eqn. (24) in B. Lucini, M. Teper, U. Wenger https://arxiv.org/abs/hep-lat/0404008.  Returns the new lattice dimensions as well as the blocked lattice configuration using latticeIndex[coord] ordering.";
simpleBlock[dir0_] := simpleBlock[dir0, Table[1, {nd}]];
simpleBlock[dir0_, anchor_] :=
    Block[
        {dims = Table[If[dir != dir0, Ceiling[latticeDimensions[[dir]]/2],
                         latticeDimensions[[dir]]], {dir, nd}],
         multiplier = Table[If[dir == dir0, 1, 2], {dir, nd}],
         oddQ = Table[If[dir != dir0, Mod[latticeDimensions[[dir]], 2] == 1],
                      {dir, nd}],
         coord, bCoord,
         result},
        result = Table[
            bCoord = latticeCoordinates[k, dims];
            coord = wrapIt[(bCoord-1)*multiplier + anchor,
                            latticeDimensions];
            If[dir1 != dir0,
               (* Just use the last link if the lattice dimension is odd. *)
               If[oddQ[[dir1]] && bCoord[[dir1]] == dims[[dir1]],
                  getLink[dir1, coord],
                  getLink[dir1, coord].getLink[dir1, shift[dir1, coord]]]],
            {dir1, nd}, {k, latticeVolume[dims]}];
        {dims, result}];

applyBlock::usage = "Create a local variable with dynamic scoping, only when needed.";
    (* Example:
       Block[{f, h, pp}, f[] = 3; h[] = 4; pp := Print[f[]]; 
             applyBlock[f, f, pp]; applyBlock[f, h, pp]]
     *)
applyBlock[var_, val_, body_] :=
    If[var === val, body, Block[{var = val}, body]];
Attributes[applyBlock] = {HoldAll};
applyBlockSteps::usage = "Apply some number of smearLinks and simpleBlock operations to a lattice, returning the new lattice dimensions as well as the resulting lattice configuration using latticeIndex[coord] ordering.";
Options[applyBlockSteps] = {"blockLevels" -> 3, "anchors" -> Automatic,
                           "smearA" -> 0.30};
timeApplyBlockSteps = 0.0;
applyBlockSteps[dir0_, OptionsPattern[]] :=
    Block[
        {indexFunction = linearSiteIndex,
         dims = latticeDimensions,
         gf = gaugeField,
         t0 = SessionTime[],
         an},
        Do[applyBlock[
            linearSiteIndex, indexFunction,
            Block[
                {latticeDimensions = dims,
                 gaugeField = gf},
                an = If[OptionValue["anchors"] === Automatic ||
                        Length[OptionValue["anchors"]] < i,
                        Table[1, {nd}],
                        OptionValue["anchors"][[i]]];
                gaugeField = smearLinks[dir0, OptionValue["smearA"]];
                {dims, gf} = simpleBlock[dir0, an];
                indexFunction = latticeIndex]],
           {i, OptionValue["blockLevels"]}];
        timeApplyBlockSteps += SessionTime[] - t0;
        {dims, gf}];

Options[polyakovLoopTallies] = Join[{
    
    }, Options[applyBlockSteps]];
(* Prevent op from matching an option *)
timePolyakovLoopTallies = 0.0;
polyakovLoopTallies["blocked", op:_?AtomQ:1, opts:OptionsPattern[]] :=
    Block[
        {debug = False, dims, gf, dir0,
         t0 = SessionTime[],
         aopts = Apply[Sequence, FilterRules[
             {opts}, Options[applyBlockSteps]]],
         tallies = Association[]},
        Do[
            {dims, gf} = applyBlockSteps[dir0, aopts];
            If[debug, Print["blocked dimensions ", dims]]; 
            Do[
                If[dir1 != dir0,
                   (* sum up over loops in the dir1 direction *)
                   Block[
                       {face = dims, uu, ntl, lt,
                        gaugeField = gf,
                        linearSiteIndex = latticeIndex,
                        latticeDimensions = dims},
                       face[[dir1]] = 1;
                       Do[
                           Block[
                               {coord = latticeCoordinates[k, face]},
                               uu = IdentityMatrix[nc];
                               Do[
                                   uu = uu.getLink[dir1, shift[dir1, coord, l]],
                                   {l, dims[[dir1]]}];
                               ntl = {dir0, dir1, coord[[dir0]]};
                               lt = Lookup[tallies, Key[ntl], {0, 0}];
                               y = stringOperator[uu, op];
                               lt += {1, y};
                               If[debug,
                                  Print["  coord=", coord, " lookup ",
                                        {ntl, KeyExistsQ[tallies, ntl], lt}]];
                               tallies[ntl] = lt],
                           {k, latticeVolume[face]}]
                   ]],
                {dir1, nd}],
            {dir0, nd}];
        timePolyakovLoopTallies += SessionTime[] - t0;
        tallies];

Options[polyakovCorrelatorTallies] = Join[{
    
    }, Options[polyakovLoopTallies]];
(* Prevent op from matching an option *)
polyakovCorrelatorTallies["blocked", op:_?AtomQ:1, opts:OptionsPattern[]] :=
    Block[
        {tallies = Association[],
         loops,
         lopts = Apply[Sequence, FilterRules[
             {opts}, Options[polyakovLoopTallies]]]},
        loops = polyakovLoopTallies["blocked", op, lopts];
        Do[
            If[dir1 != dir0,
               Do[
                   Block[
                       (* use the minimum distance including
                         looping around the lattice *)
                       {r = Min[k1 - k2,
                                    k2 + latticeDimensions[[dir0]] - k1],
                        ntl, lt, y1, y2, z},
                       y1 = loops[{dir0, dir1, k1}];
                       y2 = loops[{dir0, dir1, k2}];
                       Assert[y1[[1]] == y2[[1]]];
                       ntl = {r, latticeDimensions[[dir0]],
                              latticeDimensions[[dir1]]};
                       lt = Lookup[tallies, Key[ntl], {0, 0, 0}];
                       (* The imaginary part cancels out when averaging over slices *)
                       z = Re[Conjugate[y1[[2]]]*y2[[2]]]/(y1[[1]]*y2[[1]]);
                       lt += {1, z, z^2};
                       tallies[ntl] = lt],
                   {k1, latticeDimensions[[dir0]]},
                   {k2, k1}]],
            {dir0, nd}, {dir1, nd}];
        tallies];


polyakovLoop[dir_, anchor_, op_] :=
  (* Simple loop. *)
  Block[{u = IdentityMatrix[nc], debug = False},
   If[debug, Print["polyakovLoop lattice dimensions ", latticeDimensions]];
   Do[u = u.getLink[dir, shift[dir, anchor, i - 1]];
      If[debug, Print["   -- ",{i, Chop[stringOperator[u,1]]}]],
      {i, latticeDimensions[[dir]]}]; stringOperator[u, op]];

(* Calculate the distribution of paths across a given slice. *)
polyakovLoopNorm::usage = "See notebook \"strings.nb\" for explanation.";
Clear[polyakovLoopNorm];
polyakovLoopNorm[ll_, width_, x_] :=
    polyakovLoopNorm[ll, width, x] =
 Block[{
     ww = If[x==0,
             IdentityMatrix[width+1],
             MatrixPower[Table[x^Abs[k1 - k2],
                 {k1, width + 1}, {k2, width + 1}], ll]]},
       Simplify[Tr[ww]]];
polyakovLoop[dir0_, dir1_Integer, anchor_List, width_, expAlpha_, op_] :=
    (* Construct a strip that winds once around the lattice in direction \
      dir0 and has width w in the transverse direction dir1.  
      Calculate loop by averaging across all closed paths in the strip \
      that wind once around the lattice, weighted by
                    e^{-alpha (length of path)}.
      To make this calculationally tractable, consider only paths
      that don't move backward in the dir0 direction. *)
    (* Tests:  
      U_dir0 = const, U_dir1 = 1 should give same result as bare Polyakov loop. 
      U_dir0 = 1, U_dir1 = const should give vacuum result.
      Applying a random gauge transform should not affect results.
      Moving anchor in dir0 direction should not affect results.
      The case expAlpha=1 should be identical the old code (below). *)
 Block[{aa = anchor, slice, cc},
   cc = polyakovLoopNorm[
       latticeDimensions[[dir0]], width, expAlpha];
   (* We will apply the normalization at each dir0 step. *)
   cc = 1/cc^(1.0/latticeDimensions[[dir0]]);
   (* Choose slice containing anchor and for each dir0-
   direction link coming out of the 
   slice, look at paths containing that link, 
   averaging across links. *)
   Sum[
    slice = Table[If[k == i, IdentityMatrix[nc], zeroMatrix[]],
                  {k, width + 1}];
    (* Now, iterate through the longitudinal direction, 
    updating slice. *)
    Do[
     If[False, Print[{i,j,stringOperator[slice[[i]],op]}]];
     (* First apply the longitudinal links *)
     Do[slice[[k]] = 
        slice[[k]].getLink[dir0, shift[dir1, aa, k - 1]],
        {k, width + 1}];
     aa = shift[dir0, anchor, j];
     (* Find new slice values, multiplying by the normalization. *)
     slice = cc*(slice +
         (* Go left to right *) 
         Block[{f = zeroMatrix[]}, 
          Table[If[k > 1, 
            f = expAlpha*(f + slice[[k - 1]]).getLink[dir1, 
                  shift[dir1, aa, k - 2]]]; f,
                {k, width + 1}]] +
         (* Go right to left *)
         Block[{f = zeroMatrix[]}, 
          Reverse[Table[
            If[k < width + 1, 
             f = expAlpha*(f + slice[[k + 1]]).ConjugateTranspose[
                 getLink[dir1, shift[dir1, aa, k - 1]]]]; 
            f, {k, width + 1, 1, -1}]]]),
        {j, latticeDimensions[[dir0]]}];
    (* Finally, evaluate at starting point *)
    stringOperator[slice[[i]], op], {i, width + 1}]
 ]/;expAlpha =!=None;

polyakovLoop[dir0_, dir1_, anchor_, width1_, None, op_] :=
  (* Older version.  This is equivalent to the above in the
    case expAlpha = 1.
    Smeared over a 1-dimensional strip of a given width. 
    Normalization is so that the result is the same as the simple loop
    for lattice links = 1. *)
  Block[{aa = anchor, slice, tmpSlice, uu},
   slice = Table[
     If[i == 1, IdentityMatrix[nc], zeroMatrix[]], {i, width1 + 1}];
   Do[aa[[dir0]] = k; tmpSlice = slice;
    uu = zeroMatrix[];
    Do[uu += slice[[i - 1]];
     uu = uu.getLink[dir1, shift[dir1, aa, i - 2]];
     tmpSlice[[i]] += uu, {i, 2, width1 + 1}]; uu = zeroMatrix[];
    Do[uu += slice[[i + 1]];
     uu = uu.ConjugateTranspose[
        getLink[dir1, shift[dir1, aa, i - 1]]];
     tmpSlice[[i]] += uu, {i, width1, 1, -1}];
    (* Normalize to match simple Polyakov loop normalization *)   
    tmpSlice *= 1/(width1 + 1);
    If[False, Print[{k, Map[Tr, slice], Map[Tr, tmpSlice]}]];
    (* Move strip across links in the dir0 direction. *)
    Do[
     slice[[i]] =
     tmpSlice[[i]].getLink[dir0, shift[dir1, aa, i - 1]],
     {i, width1 + 1}], {k, latticeDimensions[[dir0]]}];
   aa[[dir0]] = 1;
   uu = slice[[width1 + 1]];
   Do[uu = uu.ConjugateTranspose[
       getLink[dir1, shift[dir1, aa, i - 1]]];
      uu += slice[[i]], {i, width1, 1, -1}]; stringOperator[uu, op]];
(*  *** Unfinished *** *)
polyakovLoop[dir0_, anchor_List, width_, expAlpha_, op_] :=
    (* Construct a D-dimensional box that winds once around the lattice
       in direction dir0 and has width w in the D-1 transverse directions.  
  Calculate the loop by averaging across all closed paths in the box that \
wind once around the lattice, weighted by e^{- 
  alpha (length of path)}.  To make this calculationally tractable, 
  consider only paths that don't move backward in the dir0 direction \
and have minimal length in each transverse hyper-plane ("slice").  *)
    Block[{coords, 
    sliceDimensions = Table[If[dir == dir0, 1, width + 1], {dir, nd}],
    sliceVolume = (width + 1)^(nd - 1), aa = anchor, vec, slice, cc},
    Print["Unfinished"]; Abort[];
   (* To start with, 
     we need to calculate the distribution of paths across any slice. *)
      ww = Table[
     expAlpha^Abs[k1 - k2], {k1, width + 1}, {k2, width + 1}];
   (* Find eigenpair with positive eigenvector *)
   {val, vec} = 
    SelectFirst[Eigensystem[ww], Apply[Equal, Sign[#[[2]]]]];
   vec = vec/Norm[vec];
   (* Normalization for slice update. *)
   cc = LinearSolve[ww, Table[1, {width + 1}]];
   (* Distribute initial identity over slice *)
   slice = Block[{f = IdentityMatrix[nc]}, 
     Table[If[k > 1, f = f.getLink[dir1, shift[dir1, anchor, k - 1]]];
       f vec[[k]], {k, width + 1}]];
   (* Now, iterate through the longitudinal direction, 
   updating slice. *)
   Do[
    (* First apply the longitudinal links *)
    Do[slice[[j]] = 
      slice[[j]].getLink[dir0, shift[dir1, aa, k - 1]], {k, 1, 
      width + 1}]; aa = shift[dir0, anchor, j];(* 
    Find new slice values *)
    slice = cc*(slice +
        (* Go left to right *) 
        Block[{f = zeroMatrix[]}, 
         Table[If[k > 1, 
           f = expAlpha*(f + slice[[k - 1]]).getLink[dir1, 
               shift[dir1, aa, k - 1]]]; f, {k, width + 1}]] +
        (* Go right to left *)
        Block[{f = zeroMatrix[]}, 
         Reverse[Table[
           If[k < width + 1, 
            f = expAlpha*(f + slice[[k + 1]]).ConjugateTranspose[
                getLink[dir1, shift[dir1, aa, k - 1]]]]; 
           f, {k, width + 1, 1, -1}]]]), {j, 
     latticeDimensions[[dir0]]}];
   (* Finally, move back to anchor point *)
   Block[{f}, 
    Do[If[k < width + 1, 
      f = slice[[k]] + 
        f.ConjugateTranspose[
          getLink[dir1, shift[dir1, anchor, k - 1]]], 
      f = slice[[k]]], {k, width + 1, 1, -1}]; 
    stringOperator[f, op]]];

loopThickness[] := beta Sqrt[2 Pi]/nc^2;

polyakovLoopTallies["simple", op_:"1"] :=
 Block[{tallies = Association[]},
  Do[
   Block[{face = latticeDimensions, nf, ta,
          lt = Lookup[tallies, latticeDimensions[[dir]], {}]},
         face[[dir]] = 1;
         nf = Apply[Times, face];
         ta = ParallelTable[
             polyakovLoop[dir, latticeCoordinates[k, face], op],
             {k, nf}];
         tallies[latticeDimensions[[dir]]] = Join[lt, ta]],
   {dir, nd}];
  tallies];
polyakovLoopTallies["smeared", width1_, x_, op_:"1"] :=
 Block[{tallies = Association[]},
  Do[
   If[dir0 != dir1,
     Block[{face = latticeDimensions, skip, nf, ta,
             lt = Lookup[tallies, latticeDimensions[[dir0]], {}]},
       face[[dir0]] = 1; skip = Max[1, width1 - 1];
       nf = Apply[Times, face];
       ta = ParallelTable[
           Block[{coords = latticeCoordinates[k, face]},
                 If[Mod[coords[[dir1]] - 1, skip] == 0,
                    polyakovLoop[dir0, dir1, coords, width1, x, op],
                    Nothing]],
           {k, nf}];
       tallies[latticeDimensions[[dir0]]] = Join[lt, ta]],
      Nothing],
   {dir0, nd}, {dir1, nd}];
  tallies];


vectorOrder[x_, y_] :=
    Block[{z = NumericalOrder[Total[x^2], Total[y^2]]},
          If[z == 0, Order[Sort[Abs[x]], Sort[Abs[y]]], z]];
makeVertices[{f_}] := If[f>1, {{-f}, {0}, {f}}, {{0}}];
makeVertices[{f_, rest__}] :=
    Flatten[Map[
        If[f>1, {Prepend[#, -f], Prepend[#, 0], Prepend[#, f]},
           {Prepend[#, 0]}]&,
        makeVertices[{rest}]], 1];
aFirstCase::usage = "Return first value of an association whose key matches the pattern.";
Attributes[aFirstCase] = {HoldFirst};
aFirstCase[assoc_, pattern_] :=
  assoc[FirstCase[Keys[assoc], pattern]];
aCases::usage = "Return values of an association whose key matches the pattern.";
Attributes[aCases] = {HoldFirst};
aCases[assoc_, pattern_] :=
    Map[assoc, Cases[Keys[assoc], pattern]];

polyakovCorrelatorCount::usage = "Number of distinct correlators, taking into account lattice symmetries and images wrapping around the lattice.";
polyakovCorrelatorCount[] := polyakovCorrelatorCount[latticeDimensions];
polyakovCorrelatorCount[dims_] :=
    Block[{x = Union[Table[{dims[[i]], Sort[Tally[Drop[dims, {i, i}]]]},
                           {i, Length[dims]}]],
          volume = Binomial[Floor[First[#]/2]+#[[2]], #[[2]]]&},
          Total[Apply[Times,
                      Map[volume, Map[Last, x], {2}], {1}]]];

(* The Polyakov loop correlators are tallied by longidinal
   dimension and a list of transverse images. *)
polyakovCorrelators[dir_, op_] :=
 Block[{nearest = 2^(nd-1),
        face = latticeDimensions, pp, nf, vertices},
  face[[dir]] = 1;
  vertices = makeVertices[face];
  nf = Apply[Times, face];
  pp = ParallelTable[
      polyakovLoop[dir, latticeCoordinates[k, face], op], {k, nf}];
  Merge[ParallelTable[
   Block[{tallies = Association[]},
    Do[Block[
        {x1 = latticeCoordinates[k1, face],
         x2 = latticeCoordinates[k2, face], ntL, z, lt, dx},
        (* In arXiv:hep-lat/0107007v2, Teper fits to a cosh().
          That is, he includes the direct distance plus wrapping
          any number of times around the lattice.

          We want to include off-axis contributions, as well.
          We include all distances in a D-1 dimensional cube of sources
          surrounding the point.

          ntL is the area between the two loops for various
          lattice wrappings plus the length of the loop
          (for the leading correction). *)
        dx = Mod[x1 - x2 + face/2, face] - face/2;
        ntL = Append[
            Map[
                (* Longitudinal component is always zero *)
                Rest[NumericalSort[Abs[#]]]&,
                Take[
                    Sort[Map[(#-dx)&, vertices], vectorOrder],
                    nearest]],
            latticeDimensions[[dir]]];
        lt = Lookup[tallies, Key[ntL], {0, 0, 0}];
        (* Imaginary part cancels out when averaging over face *)
        z = Re[Conjugate[pp[[k1]]] pp[[k2]]];
        lt += {1, z, z^2};
        tallies[ntL] = lt],
       {k1, kernel, nf, $KernelCount}, {k2, k1, nf}];
    tallies],
   {kernel, $KernelCount}], Total]];
polyakovCorrelators[dir0_, dir1_, width1_, x_, op_] :=
 Block[{nearest = 2^(nd-1),
        face = latticeDimensions, pp, skip, nf, vertices},
  face[[dir0]] = 1;
  vertices = makeVertices[face];
  skip = Max[1, width1];
  nf = Apply[Times, face];
  pp = ParallelTable[
      Block[{coords = latticeCoordinates[k, face]},
            If[Mod[coords[[dir1]] - 1, skip] == 0,
               polyakovLoop[dir0, dir1, coords, width1, x, op], Null]],
      {k, nf}];
  Merge[ParallelTable[
   Block[{tallies = Association[]},
    Do[If[pp[[k1]]=!=Null && pp[[k2]]=!=Null,
     Block[{x1 = latticeCoordinates[k1, face],
            x2 = latticeCoordinates[k2, face], dx, ntL, z, lt, sdx},
      (* See comment above. *)
      dx = Abs[Mod[x1 - x2 + face/2, face] - face/2];
      (* dir0 component of dx is always zero *)
      sdx = Rest[NumericalSort[dx]];
      (* Demand that angle is less than 45 degrees *)
      If[k1 == k2 || (dx[[dir1]]==First[sdx] && First[sdx]<sdx[[2]]),
          (* ntL is the area between the two loops for various
            lattice wrappings plus the length of the loop
            (for the leading correction). *)
        ntL = Append[
            Map[
                (* Longitudinal component is always zero *)
                Rest[NumericalSort[Abs[#]]]&, 
                Take[Sort[Map[(#-dx)&, vertices],
                          vectorOrder],
                     nearest]],
            latticeDimensions[[dir0]]];
         lt = Lookup[tallies, Key[ntL], {0, 0, 0}];
         (* The imaginary part cancels out when averaging over the face. *)
         z = Re[Conjugate[pp[[k1]]] pp[[k2]]];
         lt += {1, z, z^2};
         tallies[ntL] = lt]]],
       {k1, kernel, nf, $KernelCount}, {k2, k1, nf}];
    tallies],
   {kernel, $KernelCount}], Total]];

polyakovCorrelatorTallies["simple", op_:"1"] :=
    Merge[Table[polyakovCorrelators[dir, op], {dir, nd}],
          Total];
polyakovCorrelatorTallies["smeared", width_, x_, op_:"1"] :=
  Merge[Flatten[
      Table[If[dir0 != dir1,
               polyakovCorrelators[dir0, dir1, width, x, op],
            Nothing],
	    {dir0, nd}, {dir1, nd}]], Total];

talliesToAverageErrors = valueError[
    #[[2]]/#[[1]],
    (* Use sample standard deviation (with Bessel's correction). *)
    Sqrt[(Re[#[[3]]] - Re[#[[2]]]^2/#[[1]])/(#[[1]]*(#[[1]]-1))] +
    I Sqrt[(Im[#[[3]]] - Im[#[[2]]]^2/#[[1]])/(#[[1]]*(#[[1]]-1))]]&;
rescaleCorrelators[tallies_] :=
    Association[Map[
        (#->tallies[#]/aFirstCase[tallies, {{0..}, __, Last[#]}])&,
        Keys[tallies]]];

getSegment[dir_, k_Integer, l_] :=
    If[Head[gaugeSegments] === Association,
       If[KeyExistsQ[gaugeSegments, {dir, k, l}],
          gaugeSegments[{dir, k, l}],
          gaugeSegments[{dir, k, l}] = calculateSegment[dir, k, l]],
       calculateSegment[dir, k, l]];
calculateSegment[dir_, k_Integer, l_] :=
    Block[{x = latticeCoordinates[k], u = IdentityMatrix[nc]},
          Do[u = u.getLink[dir, x]; x = shift[dir, x], {l}];
          (* Inverse of latticeCoordinates *)
          {u, latticeIndex[x]}];
wilsonLoop::usage = "Planar Wilson loop.  Cache segments in gaugeSegments \
if gaugeSegments is an association.";
wilsonLoop::error = "Wilson loop error";
wilsonLoop[dir1_, dir2_, k_, l1_, l2_, op_String] :=
    Block[{u1, u2, u3, u4, x, y},
          {u1, x} = getSegment[dir1, k, l1];
          {u2, x} = getSegment[dir2, x, l2];
          {u4, y} = getSegment[dir2, k, l2];
          {u3, y} = getSegment[dir1, y, l1];
          If[x == y,
             stringOperator[u1.u2.ConjugateTranspose[u4.u3], op],
             Message[wilsonLoop::error]; $Failed]];
wilsonLoopOld::usage = "Direct calculation with no caching.
  Should give the same result as wilsonLoop.";
wilsonLoopOld[dir1_, dir2_, k_, l1_, l2_, op_String] :=
    Block[{u = IdentityMatrix[nc], uu = IdentityMatrix[nc],
           x = latticeCoordinates[k], y},
          y = x;
          Do[u = u.getLink[dir1, x]; x = shift[dir1, x], {l1}];
          Do[u = u.getLink[dir2, x]; x = shift[dir2, x], {l2}];
          Do[uu = uu.getLink[dir2, y]; y = shift[dir2, y], {l2}];
          Do[uu = uu.getLink[dir1, y]; y = shift[dir1, y], {l1}];
          If[x == y, stringOperator[u.ConjugateTranspose[uu], op],
             Message[wilsonLoop::error]; $Failed]];

wilsonLoopDistribution::usage = "Return the distribution of loop values (complex) for a given size Wilson loop.  Default is to show values where the imaginary part is > 0.  The result is indexed by the lattice dimensions of the plane enclosing the loop.";
wilsonLoopDistribution[w1_, w2_, op_String:"1", all_:False] :=
 (* Opposite loop orientations would give the complex conjugate *)
 Block[{absIm = If[op == "phases" || all,
                   #, If[Im[#] > 0, #, Conjugate[#]]]&},
  Merge[Flatten[Table[
      If[w1 < latticeDimensions[[dir1]] && w2 < latticeDimensions[[dir2]] &&
         If[w1 == w2, dir1 > dir2, dir1 != dir2],
         Association[{{w1, w2, latticeDimensions[[dir1]],
                       latticeDimensions[[dir2]]} ->
           Table[
               absIm[wilsonLoop[dir1, dir2, k, w1, w2, op]],
	       {k, latticeVolume[]}]}],
         Nothing],
      {dir1, nd}, {dir2, nd}]], Flatten[#, 2]&]];
wilsonLoopTallies::usage = "Return tallies of the Wilson loop values for a given size Wilson loop.  The result is indexed by the loop dimensions and the lattice dimensions of the plane enclosing the loop.";
wilsonLoopTallies[w1_, w2_, op_String:"1"] :=
  Merge[Flatten[Table[
      If[w1 < latticeDimensions[[dir1]] && w2 < latticeDimensions[[dir2]] &&
         If[w1 == w2, dir1 > dir2, dir1 != dir2],
         Association[{{w1, w2, latticeDimensions[[dir1]],
                       latticeDimensions[[dir2]]} ->
           Sum[
               Block[{z = Re[wilsonLoop[dir1, dir2, k, w1, w2, op]]},
                     {1, z, z^2}],
	       {k, latticeVolume[]}]}],
         Nothing],
      {dir1, nd}, {dir2, nd}]], Total];
averageWilsonLoop::usage = "Return average value for a given size Wilson loop.";
averageWilsonLoop[args__] :=
 Map[
     Block[{dist=Re[#]},
   valueError[Mean[dist], StandardDeviation[dist]/Sqrt[Length[dist]]]]&,
       wilsonLoopDistribution[args]];


sumStaples::usage = "Returns a general matrix in color space.";
sumStaples[dir1_, coords_] := Block[{total = zeroMatrix[]},
 Do[If[dir2 != dir1,
    (* forward staple *)
    total += getLink[dir2, shift[dir1, coords]].ConjugateTranspose[
       getLink[dir1, shift[dir2, coords]]].ConjugateTranspose[
	   getLink[dir2, coords]];
    (* backward staple.  Note that we simply demand that Re[ Tr[ F U]]],
    where F is the staple, gives the plaquette used in the action. 
    Thus, we are free to take the complex conjugate (charge congugate). *)
    total += ConjugateTranspose[
       getLink[dir2,
        shift[dir2, shift[dir1, coords], -1]]].ConjugateTranspose[
       getLink[dir1, shift[dir2, coords, -1]]].getLink[dir2,
	 shift[dir2, coords, -1]]], {dir2, nd}]; total];

stapleTest::usage = "Find the average plaquette using the staples
times the link.  This should be equal to the normal average plaquette.
Each plaquette contributes to 4 staples.";
stapleTest[] :=
  Chop[Sum[
     Tr[sumStaples[dir, latticeCoordinates[k]].getLink[dir,
	latticeCoordinates[k]]], {k, latticeVolume[]},
     {dir, nd}]*2/(latticeVolume[]*nd*(nd - 1)*nc*4)];

lineLinks[dir_, anchor_List] :=
  Block[{coords = anchor},
   ColumnForm[
    Table[coords[[dir]] = i;
	  getLink[dir, coords], {i, latticeDimensions[[dir]]}]]];

(*
  Two-point correlations of the plaquette operator.

  For ThinkPad, 6 cores:
      8^3, nc=3 lattice:  2.0 seconds
      16^3, nc=3 lattice: 150 seconds
 *)
Options[plaquetteCorrelationTallies] = {
    (* Charge conjugation/parity -1, 1 *)
    "chargeConjugation" -> 1};
plaquetteCorrelationTallies::usage = "Returns tallies as a function of (2*distance)^2.";
plaquetteCorrelationTallies[OptionsPattern[]] :=
 Block[{center2 = Table[If[i==#1 || i==#2, 1, 0],{i, nd}]&, pp,
   distance2 = With[{ld = latticeDimensions},
                   Compile[{{xa, _Integer, 1}, {xb, _Integer, 1}},
                           Dot[#, #]&[Mod[xa - xb + ld, 2*ld] - ld]]],
   pabList,
   orient = Compile[{{dir1a, _Integer}, {dir2a, _Integer}, {xa, _Integer, 1},
                     {dir1b, _Integer}, {dir2b, _Integer}, {xb, _Integer, 1}},
                    If[
                        dir1a == dir1b && dir2a === dir2b,
                        Which[
                            #[[dir1a]]^2 + #[[dir2a]]^2 == #.#, 1,
                            #[[dir1a]]==0 && #[[dir2a]]==0, 2,
                            True, 3]&[xa-xb],
                        4]]},
  pp = Flatten[Table[
      Block[{coords = latticeCoordinates[k]},
            {dir1, dir2, 2*coords + center2[dir1, dir2],
             If[OptionValue["chargeConjugation"] > 0, Re[#], Im[#]]&[
                 plaquette[dir1, dir2, coords]]/nc}],
      {dir1, nd-1}, {dir2, dir1+1, nd}, {k, latticeVolume[]}], 2];
  pabList = With[
      {pavg = Mean[Map[Last, pp]]},
      Compile[{{lt, _Real, 1}, {pa, _Real}, {pb, _Real}},
              Block[{corr = (pa - pavg)*(pb - pavg)},
                    lt + {1, corr, corr^2}]]];
 Merge[ParallelTable[
  Block[{tallies = Association[], sameDirQ,
         dir1a, dir1b, dir2a, dir2b, xa, xb, pa, pb, lt, dx2},
   Do[
       {dir1a, dir2a, xa, pa} = pp[[i]];
       {dir1b, dir2b, xb, pb} = pp[[j]];
       sameDirQ = orient[dir1a, dir2a, xa,
                         dir1b, dir2b, xb];
       (* Nearest distance, including wrap-around.
         Use twice the distance squared, so everything is an integer. *)
       dx2 = distance2[xa, xb];
       lt = Lookup[tallies, Key[{sameDirQ, dx2}], Table[0.0, {3}]];
       tallies[{sameDirQ, dx2}] = pabList[lt, pa, pb],
       {i, kernel, Length[pp]-1, $KernelCount},
       {j, i, Length[pp]}];
   tallies],
  {kernel, $KernelCount}], Total]];
(* From https://en.wikipedia.org/wiki/Correlation_and_dependence
  There are more correct formulas for the standard error. *)
sampleCorrelation[{q_, d_} -> x_] :=
    Block[{r = x[[2]]/x[[1]]},
          {Sqrt[d]/2, r, Sqrt[x[[3]] - x[[1]] r^2]/x[[1]], q}];
makePlaquetteCorrelations[opts:OptionsPattern[]] :=
    Map[sampleCorrelation,
        Sort[Normal[plaquetteCorrelationTallies[opts]],
             NumericalOrder[#1[[1, 2]], #2[[1, 2]]]&]];
