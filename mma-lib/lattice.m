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

plaquette[dir1_, dir2_, coords_List] := Re[Tr[
    getLink[dir1, coords].
    getLink[dir2, shift[dir1, coords]].
    ConjugateTranspose[getLink[dir1, shift[dir2, coords]]].
    ConjugateTranspose[getLink[dir2, coords]]]];
averagePlaquette[] :=
    Sum[Re[plaquette[dir1, dir2, latticeCoordinates[k]]],
	     {k, latticeVolume[]}, {dir1, 2, nd},
	     {dir2, dir1 - 1}]*2/(latticeVolume[]*nc*nd*(nd - 1));
makeRootLattice[] := makeRootLattice[gaugeField];
makeRootLattice[gf_] := ParallelMap[SUPower[#, 0.5]&, gf, {2}];
latticeDistance::usage = "Distance between two lattice configurations
 using the Euclidean metric in the tangent space of
 the first lattice, divided by the number of links.";
latticeDistance[root1_, lattice2_]:=
    Block[{y = Flatten[MapThread[
        Norm[Map[Tr, 2 Im[SUGenerators[].SULog[
            Block[{x=ConjugateTranspose[#1]}, x.#2.x]]]]]&,
            {root1, lattice2}, 2]]},
          Norm[y]/Length[y]];
(* Return a function that evaluates distance from
  lattice1.  So
       latticeDistanceFrom[lattice1][lattice2] ==
           latticeDistance[makeRootLattice[lattice1], lattice2]
  This might be a memory hog, since the module variable
  persists.  *)
latticeDistanceFrom[lattice1_]:=
    Module[{root1 = makeRootLattice[lattice1]},
           latticeDistance[root1, #]&];

makeTrivialLattice::usage =
  "All links identity for a given nd,nc,latticeDimensions.
Can add small random pertubation.";
Options[makeTrivialLattice] = {randomGauge -> False,
  randomPerturbation -> 0};
makeTrivialLattice[
  OptionsPattern[]] := (gaugeField =
   Table[If[OptionValue[randomPerturbation] > 0,
     MatrixExp[
      I Table[RandomReal[
          OptionValue[randomPerturbation]],
	      {nc^2 - 1}] .SUGenerators[]],
     IdentityMatrix[nc]], {nd}, {latticeVolume[]}];
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


polyakovLoop[dir_, anchor_] :=
  (* Simple loop. *)
  Block[{u = IdentityMatrix[nc]},
   Do[If[False, Print[{i, Tr[u]}]];
      u = u.getLink[dir, shift[dir, anchor, i - 1]],
      {i, latticeDimensions[[dir]]}]; Tr[u]];
polyakovLoop[dir0_, anchor_, {dir1_}, {width1_}] :=
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
    uu += slice[[i]], {i, width1 - 1, 1, -1}]; Tr[uu]];

loopThickness[] := beta Sqrt[2 Pi]/nc^2;

(* The Polyakov loop correlators are tallied as
  a function of the area enclosed by the two loops (ntl).
  This allows easy generalization to non-cubic lattices. *)
SetAttributes[polyakovLoopAdd, HoldFirst];
polyakovLoopAdd[tallies_, dir_] :=
 Block[{face = latticeDimensions, pp, nf}, face[[dir]] = 1;
  nf = Apply[Times, face];
  pp = Table[polyakovLoop[dir, latticeCoordinates[k, face]], {k, nf}];
   Do[Block[{x1 = latticeCoordinates[k1, face],
             x2 = latticeCoordinates[k2, face], ntL, z, lt},
    (* Here, we use the shortest distance, including wrapping
      around the lattice.  In arXiv:hep-lat/0107007v2, Teper
      fits to a cosh().  That is, he includes the direct distance
      plus wrapping any number of times around the lattice.
      However, this makes things complicated if we want to
      include off-axis separation or contributions from wrappings in
      the other direction transverse to the Polyakov loop.

      ntL is the area enclosed by the Polyakov loops. *)
    ntL = latticeDimensions[[dir]]*
          Norm[Mod[x1 - x2 + face/2, face] - face/2];
    lt = Lookup[tallies, ntL, {0, 0, 0}]; lt[[1]] += 1;
    (* Imaginary part cancels out when averaging over face *)
    z = Re[pp[[k1]] Conjugate[pp[[k2]]]]; lt[[2]] += z;
    lt[[3]] += z^2; tallies[ntL] = lt], {k1, nf - 1},
      {k2, k1 + 1, nf}]];
polyakovLoopAdd[tallies_, dir0_, dir1_] :=
 Block[{face = latticeDimensions, pp,
   width1 = Floor[loopThickness[] + 0.5] + 1, skip, nf},
  face[[dir0]] = 1; skip = Max[1, width1 - 1];
  nf = Apply[Times, face];
  pp = Table[
    Block[{coords = latticeCoordinates[k, face]},
     If[Mod[coords[[dir1]] - 1, skip] == 0,
      polyakovLoop[dir0, coords, {dir1}, {width1}], Null]], {k, nf}];
  Do[If[NumberQ[pp[[k1]]] && NumberQ[pp[[k2]]],
    Block[{x1 = latticeCoordinates[k1, face],
           x2 = latticeCoordinates[k2, face], dx, ntL, z, lt},
     (* See comment above. *)
     dx = Mod[x1 - x2 + face/2, face] - face/2;
     If[2*dx[[dir1]]^2 < dx.dx,
        (* ntL is the area enclosed by the Polyakov loops. *)
        ntL = latticeDimensions[[dir0]]*Norm[dx];
        lt = Lookup[tallies, ntL, {0, 0, 0}]; lt[[1]] += 1;
      (* The imaginary part cancels out when averaging over the face. *)
      z = Re[pp[[k1]] Conjugate[pp[[k2]]]]; lt[[2]] += z;
      lt[[3]] += z^2; tallies[ntL] = lt]]],
     {k1, nf - 1}, {k2, k1 + 1, nf}]];

polyakovLoopTallies["simple"] :=
  Block[{tallies = Association[{}]},
   Do[polyakovLoopAdd[tallies, dir], {dir, nd}]; tallies];
polyakovLoopTallies["smeared"] :=
  Block[{tallies = Association[{}]},
	Do[If[dir0 != dir1, polyakovLoopAdd[tallies, dir0, dir1]],
	   {dir0, nd}, {dir1, nd}]; tallies];

talliesToAverageErrors[tallies_] :=
 Map[{#[[2]]/#[[1]],
    Sqrt[(#[[1]] #[[3]] - #[[2]]^2)/#[[1]]]/#[[1]]}&, tallies];

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


setAxialGauge::usage = "Set gauge for the current lattice configuration.
This should be exact, so the average plaquette should be unchanged.";
setAxialGauge[dir_] :=
 Do[Block[{coords = latticeCoordinates[i], debug = False},
   If[coords[[dir]] == 1,
    Block[{v = IdentityMatrix[nc], vt, ct},
     (* Initial values *)
     If[False,
	Do[Print[{j, getLink[dir, shift[dir, coords, j - 1]]}],
	   {j, latticeDimensions[[dir]]}]];
     Do[v = v.getLink[dir, shift[dir, coords, j - 1]],
	{j, latticeDimensions[[dir]]}];
     If[debug, Print["Full product:", v]];
     v = SUPower[v, 1/latticeDimensions[[dir]]];
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
plotActionLimits[] := {0, Min[2 n^2/3, 15]};
plotActionPlane[dir1_, dir2_, anchor_, range_] :=
 Block[{coords = anchor},
  ListPlot3D[
   Flatten[Table[coords[[dir1]] = x1;
     coords[[dir2]] = x2; {x1 + 1/2, x2 + 1/2,
	beta (1 - plaquette[dir1, dir2, wrapIt[coords]]/nc)},
		 {x1, 0, latticeDimensions[[dir1]]},
		 {x2, 0, latticeDimensions[[dir2]]}],
	   1], PlotRange -> {{1/2, latticeDimensions[[dir1]] + 1/2},
			     {1/2, latticeDimensions[[dir2]] + 1/2},
			     range}]];
plotActionSides[dir0_, dir1_, dir2_, anchor_, range_] :=
 Block[{coords = anchor},
  ListPlot3D[
   Flatten[Table[coords[[dir1]] = x1;
     coords[[dir2]] = x2; {x1 + If[i == 1, 1/2, 0],
      x2 + If[i == 1, 0, 1/2],
      beta (1 -
         plaquette[dir0, If[i == 1, dir1, dir2], wrapIt[coords]]/nc)},
		 {x1, 0, latticeDimensions[[dir1]]},
		 {x2, 0, latticeDimensions[[dir2]]}, {i, 2}], 2],
   PlotRange -> {{1/2, latticeDimensions[[dir1]] + 1/2}, {1/2,
      latticeDimensions[[dir2]] + 1/2}, range}]];
plotAction::usage = "Animate the action density of a two dimensional slices of the lattice in directions {dir1, dir2} as a function of dir0.  One can also specify an anchor point for the remaining nd-3 dimensions.";
Options[plotAction] = {anchor :> Table[1, {nd}], PlotRange :> {0, 2 nc^2/3 }};
plotAction[dir0_, dir1_, dir2_, OptionsPattern[]] :=
 Block[{coords = OptionValue[anchor]},
       ListAnimate[
           Flatten[Table[coords[[dir0]] = x0;
               {plotActionPlane[dir1, dir2, coords, OptionValue[PlotRange]],
	        plotActionSides[dir0, dir1, dir2, coords,
                                OptionValue[PlotRange]]},
		         {x0, latticeDimensions[[dir0]]}], 1]]]

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
plaquetteCorrelationTallies::usage = "Returns tallies as a function of (2*distance)^2.";
plaquetteCorrelationTallies[] :=
 Block[{center2 = Table[If[i==#1 || i==#2, 1, 0],{i, nd}]&, pp,
   distance2 = With[{ld = latticeDimensions},
                   Compile[{{xa, _Integer, 1}, {xb, _Integer, 1}},
                           Dot[#, #]&[Mod[xa - xb + ld, 2*ld] - ld]]],
   pabList = Compile[{{lt, _Real, 1}, {pa, _Real}, {pb, _Real}},
                     lt + {1, pa, pb, pa^2, pb^2, pa*pb}]},
  pp = Flatten[Table[
      Block[{coords = latticeCoordinates[k]},
            {2*coords + center2[dir1, dir2],
             plaquette[dir1, dir2, coords]/nc}],
      {dir1, nd-1}, {dir2, dir1+1, nd}, {k, latticeVolume[]}], 2];
 Merge[ParallelTable[
  Block[{tallies = Association[{}], xa, xb, pa, pb, lt, dx2},
   Do[
       {xa, pa} = pp[[i]];
       {xb, pb} = pp[[j]];
       (* Nearest distance, including wrap-around.
         Use twice the distance squared, so everything is an integer. *)
       dx2 = distance2[xa, xb];
       lt = Lookup[tallies, dx2, Table[0.0, {6}]];
       tallies[dx2] = pabList[lt, pa, pb],
       {i, kernel, Length[pp]-1, $KernelCount},
       {j, i+1, Length[pp]}];
   tallies],
  {kernel, $KernelCount}], Total]];
(* From https://en.wikipedia.org/wiki/Correlation_and_dependence
  There are more correct formulas for the standard error. *)
sampleCorrelation[d_ -> x_] :=
    Block[{r = (x[[1]]*x[[6]] - x[[2]]*x[[3]])/
               Sqrt[(x[[1]]*x[[4]] - x[[2]]^2) (x[[1]]*x[[5]] - x[[3]]^2)]},
          {Sqrt[d]/2, r, Sqrt[(1-r^2)/(x[[1]]-2)]}];
plaquetteCorrelations[] :=
    Map[sampleCorrelation,
        Sort[Normal[plaquetteCorrelationTallies[]],
             Order[N[First[#1]], N[First[#2]]] &]];
plotPlaquetteCorrelations[corr_, opts:OptionsPattern[]] := 
   (* Bug in ErrorListPlot:  x-values must be numerical *)
 ErrorListPlot[Map[{{N[#[[1]]], #[[2]]}, ErrorBar[ #[[3]]]} &, corr], 
  PlotRange -> {{0, All}, All}, Axes -> False, Frame -> True, 
  FrameLabel -> {"distance", "correlation"}, 
  Epilog -> {Dashing[0.01], 
             Line[{{0, 0}, {1/2 + Max[Map[First, corr]], 0}}]}, opts];
