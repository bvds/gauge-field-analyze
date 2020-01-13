(*
      Quantities associated with a Wilson action lattice

   A lattice configuration is defined via the global variables:
     nc:  number of colors
     nd:  number of dimensions
     linearSiteIndex:  Map of lattice coordinates onto integers
     latticeDimensions:  Array containing the size of the lattice
     gaugeField:  Array dimensions {nd, latticeVolume[], nc, nc} containing
                  the gauge fields
*)

getLink[dir_, coords_List] :=
  gaugeField[[dir, linearSiteIndex[coords]]];
getLink[rgf_][dir_, coords_] := rgf[[dir, linearSiteIndex[coords]]];
setLink[dir_, coords_List, mat_] :=
 gaugeField[[dir, linearSiteIndex[coords]]] = mat;
latticeVolume[] := latticeVolume[latticeDimensions];
latticeVolume[dimensions_] := Apply[Times, dimensions];

latticeCoordinates::usage = "Given an integer-valued site index,
return lattice coordinates.  This is not necessarily the inverse of
linearSiteIndex.  It is simply a way to iterate over the whole lattice.";
latticeCoordinates[i_] := latticeCoordinates[i, latticeDimensions]
latticeCoordinates[ii_, dimensions_] :=
 Block[{i = ii - 1, j},
  Map[({i, j} = QuotientRemainder[i, #]; j + 1)&, dimensions]];

latticeIndex::usage = "Inverse of latticeCoordinates.";
latticeIndex[coords_List] := latticeIndex[coords, latticeDimensions];
latticeIndex[coords_List, dimensions_] :=
 Block[{i = 1, delta = 1},
  Do[i += (coords[[k]] - 1)*delta;
   delta *= dimensions[[k]], {k, Length[dimensions]}]; i];

shift::usage = "Wrap-around, assuming periodic boundary conditions.";
shift[dir_, coordIn_List, size_:1] :=
  Block[{cc = coordIn},
	cc[[dir]] = 1 + Mod[cc[[dir]] + size - 1, latticeDimensions[[dir]]];
	cc];

plaquette[dir1_, dir2_, coords_List] := Tr[
    getLink[dir1, coords].
    getLink[dir2, shift[dir1, coords]].
    ConjugateTranspose[getLink[dir1, shift[dir2, coords]]].
    ConjugateTranspose[getLink[dir2, coords]]];
averagePlaquette[] :=
    Sum[Re[plaquette[dir1, dir2, latticeCoordinates[k]]],
	     {k, latticeVolume[]}, {dir1, 2, nd},
	     {dir2, dir1 - 1}]*2/(latticeVolume[]*nc*nd*(nd - 1));
makeRootLattice[] := makeRootLattice[gaugeField];
makeRootLattice[gf_] := ParallelMap[SUPower[#, 0.5]&, gf, {2}];
latticeDistance::usage = "Distance between two lattice configurations
 using the Euclidean metric in the tangent space of
 the first lattice, divided by the number of links.";
latticeDistance[lattice1_, lattice2_]:=
    Block[{y = Flatten[MapThread[
        (* These are all equivalent *)
        Which[
            True,
            First[SUNorm[LinearSolve[##]]],
            False,
            Sqrt[2] Norm[Flatten[SULog[LinearSolve[##]]]],
            False,
            Sqrt[2 Tr[#.ConjugateTranspose[#]]]&[SULog[LinearSolve[##]]],
            False,
            (* The Im[] is to remove any floating point errors *)
            2 Norm[Map[Tr, Im[SUGenerators[].SULog[LinearSolve[##]]]]]]&,
            {lattice1, lattice2}, 2]]},
          Norm[y]/Sqrt[Length[y]]];
latticeNorm::usage = "SUNorm[] averaged over the lattice links.";
latticeNorm[opts:OptionsPattern[]] := latticeNorm[gaugeField, opts];
latticeNorm[lattice_, opts:OptionsPattern[]]:=
    Block[{y = Flatten[Map[
        First[SUNorm[#, opts]]&,
             lattice, {2}]]},
          Norm[y]/Sqrt[Length[y]]];

makeTrivialLattice::usage =
  "All links identity for a given nd,nc,latticeDimensions.
Can add small random pertubation or choose and element of the center.";
Options[makeTrivialLattice] = {randomGauge -> False,
  randomCenter -> False,
  randomPerturbation -> 0};
makeTrivialLattice[
  OptionsPattern[]] := (gaugeField =
   Table[If[OptionValue[randomPerturbation] > 0,
     MatrixExp[
      I Table[RandomReal[
          OptionValue[randomPerturbation]],
	      {nc^2 - 1}] .SUGenerators[]],
            IdentityMatrix[nc]]*
         If[OptionValue[randomCenter],
            Exp[I 2 Pi RandomInteger[nc-1]/nc], 1],
         {nd}, {latticeVolume[]}];
  Clear[linearSiteIndex];
  Do[linearSiteIndex[latticeCoordinates[i]] = i,
     {i, latticeVolume[]}];
  If[OptionValue[randomGauge],
   Do[Block[{coords = latticeCoordinates[i], u = randomSUMatrix[]},
     Do[setLink[dir, coords, u.getLink[dir, coords]];
      setLink[dir, shift[dir, coords, -1],
	      getLink[dir, shift[dir, coords, -1]].
		     ConjugateTranspose[u]], {dir, nd}]],
      {i, latticeVolume[]}]]);

polyakovLoop[dir_, anchor_, op_] :=
  (* Simple loop. *)
  Block[{u = IdentityMatrix[nc]},
   Do[If[False, Print[{i, Tr[u]}]];
      u = u.getLink[dir, shift[dir, anchor, i - 1]],
      {i, latticeDimensions[[dir]]}]; stringOperator[u, op]];
polyakovLoop[dir0_, anchor_, {dir1_}, {width1_}, op_] :=
  (* Smeared over a 1-dimensional strip of a given width. 
  Normalization is so that the result is the same as the simple loop
  for lattice links = 1. *)
  Block[{aa = anchor, slice, tmpSlice, uu},
   slice = Table[
     If[i == 1, IdentityMatrix[nc], zeroMatrix[]], {i, width1}];
   Do[aa[[dir0]] = k; tmpSlice = slice;
    uu = zeroMatrix[];
    Do[uu += slice[[i - 1]];
     uu = uu.getLink[dir1, shift[dir1, aa, i - 2]];
     tmpSlice[[i]] += uu, {i, 2, width1}]; uu = zeroMatrix[];
    Do[uu += slice[[i + 1]];
     uu = uu.ConjugateTranspose[
        getLink[dir1, shift[dir1, aa, i - 1]]];
     tmpSlice[[i]] += uu, {i, width1 - 1, 1, -1}];
    (* Normalize to match simple Polyakov loop normalization *)   
    tmpSlice *= 1/width1;
    If[False, Print[{k, Map[Tr, slice], Map[Tr, tmpSlice]}]];
    (* Move strip across links in the dir0 direction. *)
    Do[
     slice[[i]] =
     tmpSlice[[i]].getLink[dir0, shift[dir1, aa, i - 1]],
     {i, width1}], {k, latticeDimensions[[dir0]]}];
   aa[[dir0]] = 1;
   uu = slice[[width1]];
   Do[uu = uu.ConjugateTranspose[
       getLink[dir1, shift[dir1, aa, i - 1]]];
    uu += slice[[i]], {i, width1 - 1, 1, -1}]; stringOperator[uu, op]];

loopThickness[] := beta Sqrt[2 Pi]/nc^2;

(* The Polyakov loop correlators are tallied as
  a function of the area enclosed by the two loops (ntL).
  This allows easy generalization to non-cubic lattices. *)
SetAttributes[polyakovLoopAdd, HoldFirst];
polyakovLoopAdd[tallies_, dir_, op_] :=
 Block[{face = latticeDimensions, pp, nf}, face[[dir]] = 1;
  nf = Apply[Times, face];
  pp = Table[polyakovLoop[dir, latticeCoordinates[k, face], op], {k, nf}];
  Do[Block[
      {x1 = latticeCoordinates[k1, face],
       x2 = latticeCoordinates[k2, face], ntL, z, lt},
      (* Here, we use the shortest distance, including wrapping
        around the lattice.  In arXiv:hep-lat/0107007v2, Teper
        fits to a cosh().  That is, he includes the direct distance
        plus wrapping any number of times around the lattice.
        However, this makes things complicated if we want to
        include off-axis separation or contributions from wrappings in
        the other direction transverse to the Polyakov loop.

        ntL is the area enclosed by the Polyakov loop. *)
      ntL = latticeDimensions[[dir]]*
            Norm[Mod[x1 - x2 + face/2, face] - face/2];
      lt = Lookup[tallies, ntL, {0, 0, 0}];
      If[Head[op] === String,
         (* Imaginary part cancels out when averaging over face *)
         z = Re[Conjugate[pp[[k1]]] pp[[k2]]];
         lt += {1, z, z^2},
         z = Conjugate[pp[[k1, 1]]] pp[[k2, 2]];
         lt += {1, z, Re[z]^2 + I*Im[z]^2};
         z = Conjugate[pp[[k2, 1]]] pp[[k1, 2]];
         lt += {1, z, Re[z]^2 + I*Im[z]^2}];
      tallies[ntL] = lt],
     {k1, nf - 1}, {k2, k1 + 1, nf}]];
polyakovLoopAdd[tallies_, dir0_, dir1_, op_] :=
 Block[{
  face = latticeDimensions, pp,
  width1 = Floor[loopThickness[] + 0.5] + 1, skip, nf},
       face[[dir0]] = 1; skip = Max[1, width1 - 1];
       nf = Apply[Times, face];
  pp = Table[
      Block[{coords = latticeCoordinates[k, face]},
            If[Mod[coords[[dir1]] - 1, skip] == 0,
               polyakovLoop[dir0, coords, {dir1}, {width1}, op], Null]],
      {k, nf}];
  Do[If[pp[[k1]]=!=Null && pp[[k2]]=!=Null,
    Block[{x1 = latticeCoordinates[k1, face],
           x2 = latticeCoordinates[k2, face], dx, ntL, z, lt},
     (* See comment above. *)
     dx = Mod[x1 - x2 + face/2, face] - face/2;
     If[2*dx[[dir1]]^2 < dx.dx,
        (* ntL is the area enclosed by the Polyakov loop. *)
        ntL = latticeDimensions[[dir0]]*Norm[dx];
        lt = Lookup[tallies, ntL, {0, 0, 0}];
        If[Head[op]===String,
           (* The imaginary part cancels out when averaging over the face. *)
           z = Re[Conjugate[pp[[k1]]] pp[[k2]]];
           lt += {1, z, z^2},
           z = Conjugate[pp[[k1, 1]]] pp[[k2, 2]];
           lt += {1, z, Re[z]^2 + I*Im[z]^2};
           z = Conjugate[pp[[k2, 1]]] pp[[k1, 2]];
           lt += {1, z, Re[z]^2 + I*Im[z]^2}];
        tallies[ntL] = lt]]],
     {k1, nf-1}, {k2, k1 + 1, nf}]];

polyakovLoopTallies["simple", op_:"1"] :=
  Block[{tallies = Association[{}]},
   Do[polyakovLoopAdd[tallies, dir, op], {dir, nd}]; tallies];
polyakovLoopTallies["smeared", op_:"1"] :=
  Block[{tallies = Association[{}]},
	Do[If[dir0 != dir1, polyakovLoopAdd[tallies, dir0, dir1, op]],
	   {dir0, nd}, {dir1, nd}]; tallies];

talliesToAverageErrors[tallies_] :=
    Map[{#[[2]]/#[[1]],
         Sqrt[Re[#[[3]]] - Re[#[[2]]]^2/#[[1]]]/#[[1]]
         + I Sqrt[Im[#[[3]]] - Im[#[[2]]]^2/#[[1]]]/#[[1]]}&, tallies];

Options[exponentialModel] = {printResult -> False};
exponentialModel[tallyData_, OptionsPattern[]] :=
    Block[(* Protect against any global definitions of model paramters. *)
	{ff, a2sigma, norm, c, x},
  ff = NonlinearModelFit[
    Map[{#[[1]], #[[2, 1]]}&, Normal[tallyData]],
    norm Exp[-a2sigma x] + c, {{norm, 1}, {a2sigma, 0.01}, {c, 0}}, x,
    VarianceEstimatorFunction -> (1&),
    Weights -> Map[(1/#[[2, 2]]^2)&, Normal[tallyData]]];
  If[OptionValue[printResult], Print[ff["ParameterTable"]]]; ff];


wilsonLoop::usage = "Planar Wilson loop.";
wilsonLoop[dir1_, dir2_, coords_List, l1_, l2_, op_String] :=
    Block[{u = IdentityMatrix[nc], x = coords},
          Do[u = u.getLink[dir1, x]; x = shift[dir1, x], {l1}];
          Do[u = u.getLink[dir2, x]; x = shift[dir2, x], {l2}];
          Do[x = shift[dir1, x, -1];
             u = u.ConjugateTranspose[getLink[dir1, x]], {l1}];
          Do[x = shift[dir2, x, -1];
             u = u.ConjugateTranspose[getLink[dir2, x]], {l2}];
          If[x != coords, Print["Wilson loop error"]; Abort[]];
          stringOperator[u, op]]
wilsonLoopDistribution::usage = "Return distribution of loop values (complex) for a given size Wilson loop.  Just show values where the imaginary part is > 0.";
wilsonLoopDistribution[l1_, l2_, op_String:"1"] :=
    (* Opposite loop orientations would give the complex conjugate *)
    Block[{absIm = If[Im[#] > 0, #, Conjugate[#]]&},
    Flatten[ParallelTable[
        {absIm[wilsonLoop[dir1, dir2, latticeCoordinates[k], l1, l2, op]],
         If[l1 != l2,
            absIm[wilsonLoop[dir1, dir2, latticeCoordinates[k], l2, l1, op]],
            Nothing]},
	{k, latticeVolume[]}, {dir1, 2, nd},
	{dir2, dir1 - 1}]]];
averageWilsonLoop::usage = "Return average value for a given size Wilson loop.";
averageWilsonLoop[l1_, l2_, op_String:"1"] :=
    (* Opposite loop orientations would cancel the imaginary part. *)
    (ParallelSum[
        Re[wilsonLoop[dir1, dir2, latticeCoordinates[k], l1, l2, op]],
        + If[l1 != l2,
             Re[wilsonLoop[dir1, dir2, latticeCoordinates[k], l2, l1, op]],
             0],
	{k, latticeVolume[]}, {dir1, 2, nd},
	{dir2, dir1 - 1}]*2/(If[l1 == l2, 1, 2]*latticeVolume[]*nd*(nd - 1)));


setAxialGauge::usage = "Set axial gauge for the current lattice configuration.
That is, links along any Polyakov loop in direction dir are constant.
In the case \"center\" -> True,  set the gauge modulo the center of the group.
For each Polyakov loop, one link will be assigned the total center
rotation for that loop.";
Options[setAxialGauge] = Options[SUPower];
setAxialGauge[dir_, opts:OptionsPattern[]] :=
 Do[Block[{coords = latticeCoordinates[i], debug = OptionValue["debug"]},
   If[coords[[dir]] == 1,
    Block[{v = IdentityMatrix[nc], vt, ct},
     (* Initial values *)
     If[False,
	Do[Print[{j, getLink[dir, shift[dir, coords, j - 1]]}],
	   {j, latticeDimensions[[dir]]}]];
     Do[v = v.getLink[dir, shift[dir, coords, j - 1]],
	{j, latticeDimensions[[dir]]}];
     If[debug, Print["Full product:", v]];
     v = SUPower[v, 1/latticeDimensions[[dir]], opts];
     If[debug, Print["Root:", v]];
     Do[
      (* Apply gauge choice to each site *)
      ct = shift[dir, coords, j - 1];
      vt = ConjugateTranspose[getLink[dir, shift[dir, ct, -1]]].v;
      Do[setLink[dd, shift[dd, ct, -1],
	(* In order to decrease the floating point error,
        just set the axial direction link to the correct value. *)
        If[dd == dir, v, getLink[dd, shift[dd, ct, -1]].vt]];
	 setLink[dd, ct, ConjugateTranspose[vt].getLink[dd, ct]],
	   {dd, nd}], {j, 2, latticeDimensions[[dir]]}];
     (* Verify that it all worked out *)
     If[debug,
	Do[Print[{j, getLink[dir, shift[dir, coords, j - 1]]}],
	   {j, latticeDimensions[[dir]]}];
      Print["Final link diff:",
	    Chop[v - getLink[dir, shift[dir, coords, -1]]]]]]]],
    {i, latticeVolume[]}];

(*
  Using "center"->False then "center"->True works substantially better.
  For 16^3, nc=3, beta=28, "center"->False:  "damping"->1.7 works best
         For saddle point, "center"->False has same behavior.
         For saddle point, "center"->False then "center"->True:
                                             "damping"->1.8 works best
 *)
Options[setMinimumGauge] = {"center" -> False, "iterations" -> 1,
                            "damping" -> 1.7, "debug" -> False};
setMimimumGauge::usage = "Set gauge to mimium average link magnitude for the current lattice configuration.  This should be exact, so the average plaquette should be unchanged.";
setMinimumGauge[OptionsPattern[]] :=
    Table[Block[{normSum = 0.0, update,
              damping = OptionValue["damping"], logLog},
   (* Take the log of a link, accumulating statistics
      on the first run through the checkerboard. *)
   logLog = Block[{x = SULog[#1, "center"->OptionValue["center"]], y},
                  If[#2==0, y = Flatten[x];
                            normSum += 2.0*Total[Re[y]^2+Im[y]^2]]; x]&;
   (* Take result from parallel computation and update gaugeField,
      accumulating statistics. *)
   update = (Apply[setLink, Take[#, 3]];
             If[#[[4]] =!= Null, normSum += #[[4]]])&;
   (* Shared variables are super-slow.  Instead, create a table
     of updated links, then update gaugeField and statistics after
     the parallel computation is completed.
     Use checkerboard to avoid conflicts in link updates. *)
   Do[Scan[update, ParallelTable[
       Block[{coords = latticeCoordinates[i], sum, gauge,
              normSum = 0.0},
             If[Mod[Total[coords], 2] == cb,
                (* Construct gauge transform *)
                sum = Sum[
                    logLog[getLink[dir, coords], cb]
                    - logLog[getLink[dir, shift[dir, coords, -1]], cb],
                    {dir, nd}];
                gauge = MatrixExp[-damping*sum/(2*nd)];
                (* Create gauge transformed links.
                 Include statistics once. *)
                Table[
                    {{dir, coords, gauge.getLink[dir, coords],
                      If[dir==1, normSum, Null]},
                     {dir, shift[dir, coords, -1],
                      getLink[dir, shift[dir, coords, -1]].
                             ConjugateTranspose[gauge], Null}},
                    {dir, nd}],
                Nothing]],
       {i, latticeVolume[]}, Method->"CoarsestGrained"], {3}], {cb, 0, 1}];
   (* The 2-norm of the initial links,
     averaging over the number of links. *)
   Sqrt[normSum/(latticeVolume[]*nd)]],
                {OptionValue["iterations"]}];

Options[setMinimumAxialGauge] = {"center" -> False,
                            "damping" -> 1.0, "debug" -> False};
setMinimumAxialGauge[dir1_, OptionsPattern[]] :=
 Block[{damping = OptionValue["damping"],
     (* Take the log of a link. *)
     logLog = SULog[#1, "center"->OptionValue["center"]]&,
     (* Take result from parallel computation and update gaugeField. *)
     update = Apply[setLink, #]&},
   (* Shared variables are super-slow.  Instead, create a table
     of updated links, then update gaugeField and statistics after
     the parallel computation is completed.
     Use checkerboard to avoid conflicts in link updates. *)
   Do[Scan[update, ParallelTable[
     Block[{coords = latticeCoordinates[i], sum, gauge, lastU, firstU},
      If[coords[[dir1]] == 1 && Mod[Total[coords], 2] == cb,
         lastU = getLink[dir1, shift[dir1, coords, -1]];
         Table[
             (* Iterate through a Polyakov loop in direction dir1 *)
             coords[[dir1]] = j;
             (* Construct gauge transform *)
             sum = Sum[
                 If[dir1 == dir2,
                    -logLog[lastU],
                    logLog[getLink[dir2, coords]]
                    - logLog[getLink[dir2, shift[dir2, coords, -1]]]],
                 {dir2, nd}];
             gauge = MatrixExp[-damping*sum/(2*nd-1)];
             (* Create gauge transformed links.
               Include statistics once. *)
             Table[
                 If[dir1 == dir2,
                    (* Order is important here:  need to use
                      lastU before updating it. *)
                    {If[j==1,
                        firstU = lastU.ConjugateTranspose[gauge];
                        Nothing,
                        {dir2, shift[dir2, coords, -1],
                         lastU.ConjugateTranspose[gauge]}],
                     If[j==latticeDimensions[[dir1]],
                        {dir2, coords, gauge.firstU},
                        lastU = gauge.getLink[dir2, coords];
                        Nothing]},
                    {{dir2, shift[dir2, coords, -1],
                      getLink[dir2, shift[dir2, coords, -1]].
                             ConjugateTranspose[gauge]},
                     {dir2, coords, gauge.getLink[dir2, coords]}}],
                 {dir2, nd}],
             {j, latticeDimensions[[dir1]]}],
         Nothing]],
     {i, latticeVolume[]}, Method->"CoarsestGrained"], {4}], {cb, 0, 1}]];

nStrategies = 10;
applyGaugeTransforms[s_] :=
    (* New direction each time setAxialGauge[] is called. *)
  Block[{lastDir = nd, dir}, dir = (lastDir = Mod[lastDir, nd] + 1) &;
    Scan[Which[
      # == 1, setAxialGauge[dir[], "center" -> False],
      # == 2, setAxialGauge[dir[], "center" -> True],
      # == 3, setMinimumGauge["center" -> False, "damping" -> 1],
      # == 4, setMinimumGauge["center" -> True, "damping" -> 1],
      # == 5, setMinimumGauge["center" -> False, "damping" -> 1.5],
      # == 6, setMinimumGauge["center" -> True, "damping" -> 1.5],
      # == 7,
      setMinimumAxialGauge[dir[], "center" -> False, "damping" -> 1],
      # == 8,
      setMinimumAxialGauge[dir[], "center" -> True, "damping" -> 1],
      # == 9,
      setMinimumAxialGauge[dir[], "center" -> False, "damping" -> 1.5],
      # == 10,
      setMinimumAxialGauge[dir[], "center" -> True, "damping" -> 1.5]
         ] &, s]; latticeNorm["center" -> True]];

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

wrapIt[coords_List] :=
    MapThread[(1 + Mod[#1 - 1, #2])&, {coords, latticeDimensions}];

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
  Block[{tallies = Association[{}], sameDirQ,
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
       {j, i+1, Length[pp]}];
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
             Order[N[#1[[1, 2]]], N[#2[[1, 2]]]]&]];

makeFaradayLattice::usage = 
  "Construct a dual lattice containing the Faraday tensor.  This \
construction only makes sense in the context of the minimum norm \
gauge."; 
makeFaradayLattice[] := 
  Block[{lgf = (-I)*ParallelMap[SULog, gaugeField, {2}], agA1, agA2, agA1s2, 
    agA2s1, agA11, agA22, faces = {{2, 3}, {3, 1}, {1, 2}}, 
    dx = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}}, 
    dualLattice = Array[Null &, {nd, latticeVolume[]}], dir1, dir2, 
    coords}, 
  Do[{dir1, dir2} = faces[[dir0]]; 
      coords = latticeCoordinates[i];
    agA1 = getLink[lgf][dir1, coords]; 
    agA2s1 = getLink[lgf][dir2, shift[dir1, coords]]; 
    agA1s2 = getLink[lgf][dir1, shift[dir2, coords]]; 
    agA2 = getLink[lgf][dir2, coords];
    agA11 = (agA1 + agA1s2)/2; agA22 = (agA2 + agA2s1)/2; 
    coords = wrapIt[coords - dx[[dir0]]];
    dualLattice[[dir0, linearSiteIndex[coords]]] =
    (agA2s1 - agA2) - (agA1s2 - agA1) + 
      I*(agA11.agA22 - agA22.agA11), {dir0, nd}, 
         {i, latticeVolume[]}]; dualLattice];
