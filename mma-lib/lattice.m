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
            First[SUNorm[#1.ConjugateTranspose[#2]]],
            False,
            First[SUNorm[LinearSolve[##]]],
            False,
            Norm[Flatten[SULog[LinearSolve[##]]]],
            False,
            Sqrt[Tr[#.ConjugateTranspose[#]]]&[SULog[LinearSolve[##]]]]&,
            {lattice1, lattice2}, 2]]},
          Norm[y]/Sqrt[Length[y]]];
latticeNorm::usage = "SUNorm[] averaged over the lattice links.";
Options[latticeNorm] = Join[Options[SUNorm], {"tally" -> False}];
latticeNorm[opts:OptionsPattern[]] := latticeNorm[gaugeField, opts];
latticeNorm[lattice_, opts:OptionsPattern[]]:=
    Block[{y, dis,
           sopts = Apply[Sequence, FilterRules[{opts}, Options[SUNorm]]]},
          {y, dis} = Transpose[Map[
              SUNorm[#, sopts]&, Flatten[lattice, 1]]];
          If[OptionValue["tally"],
             {Norm[y]/Sqrt[Length[y]], Sort[Tally[dis]]},
             Norm[y]/Sqrt[Length[y]]]];

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
 Block[{tallies = Association[{}]},
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
 Block[{tallies = Association[{}]},
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


makeVertices[{f_}] := {{0}, If[f>1, {-f}, Nothing], If[f>1, {f}, Nothing]}
makeVertices[{f_, rest__}] :=
    Flatten[Map[
        {Prepend[#, 0], If[f>1, Prepend[#, -f], Nothing],
         If[f>1, Prepend[#, f], Nothing]}&,
        makeVertices[{rest}]], 1];
aFirstCase::usage = "Return first value of an association whose key matches the pattern.";
Attributes[aFirstCase] = {HoldFirst};
aFirstCase[assoc_, pattern_] :=
  assoc[FirstCase[Keys[assoc], pattern]];
aCases::usage = "Return values of an association whose key matches the pattern.";
Attributes[aCases] = {HoldFirst};
aCases[assoc_, pattern_] :=
    Map[assoc, Cases[Keys[assoc], pattern]];

(* The Polyakov loop correlators are tallied as
  a function of the area enclosed by the two loops (ntL).
  This allows easy generalization to non-cubic lattices. *)
polyakovCorrelators[dir_, op_] :=
 Block[{nearest = Max[5, 2^(nd-1)],
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
        dx = Abs[Mod[x1 - x2 + face/2, face] - face/2];
        ntL = Append[
            Take[Sort[Map[Norm[#-dx]&, vertices], Less], nearest],
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
 Block[{nearest = Max[5, 2^(nd-1)],
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
      sdx = Rest[Sort[dx, Less]];
      (* Demand that angle is less than 45 degrees *)
      If[k1 == k2 || (dx[[dir1]]==First[sdx] && First[sdx]<sdx[[2]]),
          (* ntL is the area between the two loops for various
            lattice wrappings plus the length of the loop
            (for the leading correction). *)
         ntL = Append[
             Take[Sort[Map[Norm[#-dx]&, vertices], Less], nearest],
             latticeDimensions[[dir0]]];
         lt = Lookup[tallies, Key[ntL], {0, 0, 0}];
         (* The imaginary part cancels out when averaging over the face. *)
         z = Re[Conjugate[pp[[k1]]] pp[[k2]]];
         lt += {1, z, z^2};
         tallies[ntL] = lt]]],
       {k1, kernel, nf, $KernelCount}, {k2, k1, nf}];
    tallies],
   {kernel, 1 $KernelCount}], Total]];

polyakovCorrelatorTallies["simple", op_:"1"] :=
    Merge[Table[polyakovCorrelators[dir, op], {dir, nd}],
          Total];
polyakovCorrelatorTallies["smeared", width_, x_, op_:"1"] :=
  Merge[Flatten[
      Table[If[dir0 != dir1,
               polyakovCorrelators[dir0, dir1, width, x, op],
            Nothing],
	    {dir0, nd}, {dir1, nd}]], Total];

talliesToAverageErrors[tallies_] :=
    Map[valueError[#[[2]]/#[[1]],
         Sqrt[Re[#[[3]]] - Re[#[[2]]]^2/#[[1]]]/#[[1]]
         + I Sqrt[Im[#[[3]]] - Im[#[[2]]]^2/#[[1]]]/#[[1]]]&, tallies];
rescaleCorrelators[tallies_] :=
    Association[Map[
        (#->tallies[#]/aFirstCase[tallies, {0, __, Last[#]}])&,
        Keys[tallies]]];

(* See file "coulomb.nb" *)
selfEnergyCoulomb[w_, eps_:1] := w*(Log[w/eps] - 1);
twoSidesCoulomb[wa_, wb_] :=
    Block[{eta = Sqrt[wa^2 + wb^2]}, 
          2 wa - 2*eta + 2*wb*Log[(eta + wb)/wa]];
wilsonCoulomb[w1_, w2_, eps_:1]:=
    2*selfEnergyCoulomb[w1, eps] + 2*selfEnergyCoulomb[w2, eps] -
    twoSidesCoulomb[w1, w2] - twoSidesCoulomb[w2, w1];
wilsonCoulomb[w1_, w2_, l1_, l2_, eps_:1]:=
    wilsonCoulomb[w1, w2, eps] -
    twoSidesCoulomb[l1 - w1, w2] - twoSidesCoulomb[l2 - w2, w1];

modelChi2::usage = "Calculate chi^2 statistic from errors and differences.";
modelChi2[model_, errors_] :=
    Total[(model["FitResiduals"]/errors)^2];

stringModel::usage = "Fit to an exponential plus a constant term, including the universal string correction.  See Andreas Athenodorou, Michael Teper, https://arxiv.org/abs/1609.03873; arXiv:1302.6257v2 [hep-th] 12 Mar 2013; http://arxiv.org/abs/0912.3339";
Options[stringModel] = {printResult -> False, "lowerCutoff"->0,
                        "state" -> 0, "constantTerm" -> False};
Clear[c0, c1, cCoulomb, eEpsilon];
Format[c0]:=Subscript["c", 0];
Format[c1]:=Subscript["c", 1];
Format[cCoulomb]:=Subscript["c", "q"];
Format[cEpsilon]:=Subscript["c", "p"];
stringModel[tallyData_, OptionsPattern[]] :=
 Block[(* Protect against any global definitions of model parameters. *)
     {ff, r, ll,
      nll = Length[Union[Map[Last, Keys[tallyData]]]],
      na = Length[First[Keys[tallyData]]] - 1,
      casimir = 8 Pi (OptionValue["state"] - (nd-2)/24),
      data = Select[Normal[tallyData],
                    #[[1,1]]*#[[1,2]]>OptionValue["lowerCutoff"]&]},
  ff = NonlinearModelFit[
      Map[Append[#[[1]], #[[2, 1]]]&, data],
      (* The constrained version of NonLinearModelFit[] uses a
        different algorithm and does not converge as well.  Instead,
        add boundary terms to the fit function.

       Tried Exp[-const*r[i]] factor in Coulomb term for 12x18x18,
       but preferred sign was const<0. *)
      c0*Sum[
          Exp[-r[i]*ll*a2sigma*Sqrt[Max[0, 1 + casimir/(a2sigma*ll^2)] +
                                    10^-10] +
              If[nll>1, cEpsilon*ll, 0] - 2*cCoulomb*ll*Log[r[i]]] +
          r[i]*If[cCoulomb<0, 100*cCoulomb, 0],
          {i, na}] + If[OptionValue["constantTerm"], c1, 0],
      {{c0, 1.0}, {a2sigma, 0.016}, {cCoulomb, 0.02},
       If[nll>1, {cEpsilon, 0.}, Nothing],
       If[OptionValue["constantTerm"], {c1, 0.}, Nothing]},
      Append[Table[r[i], {i, na}], ll],
      VarianceEstimatorFunction -> (1&),
      Weights -> Map[(1/#[[2, 2]]^2)&, data]];
  If[OptionValue[printResult],
     Print[ff["ParameterTable"]];
     Print["chi^2: ", modelChi2[ff, Map[#[[2,2]]&, data]],
           " for ",
           Length[data]-Length[ff["BestFitParameters"]], " d.o.f."]];
  ff];


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
wilsonLoopDistribution::usage = "Return the distribution of loop values (complex) for a given size Wilson loop.  Default is to show values where the imaginary part is > 0.  The result is indexed by the lattice dimensions of the plane enclosing the loop.";
wilsonLoopDistribution[l1_, l2_, op_String:"1", all_:False] :=
 (* Opposite loop orientations would give the complex conjugate *)
 Block[{absIm = If[op == "phases" || all,
                   #, If[Im[#] > 0, #, Conjugate[#]]]&},
  Merge[Flatten[Table[
      If[l1<=latticeDimensions[[dir1]] && l2<=latticeDimensions[[dir2]] &&
         If[l1 == l2, dir1 > dir2, dir1 != dir2],
         Association[{{latticeDimensions[[dir1]], latticeDimensions[[dir2]]} ->
           ParallelTable[
               absIm[wilsonLoop[dir1, dir2, latticeCoordinates[k],
                                l1, l2, op]],
	       {k, latticeVolume[]}]}],
         Nothing],
      {dir1, nd}, {dir2, nd}]], Flatten[#, 2]&]];
averageWilsonLoop::usage = "Return average value for a given size Wilson loop.";
averageWilsonLoop[args__] :=
 Map[
     Block[{dist=Re[#]},
   valueError[Mean[dist], StandardDeviation[dist]/Sqrt[Length[dist]]]]&,
       wilsonLoopDistribution[args]];

setGlobalGauge::usage = "Perform global color rotation. u is a unitary matrix.";
setGlobalGauge[u_] :=
    (gaugeField = Map[u.#.ConjugateTranspose[u]&, gaugeField, {2}]);

setRandomGauge::usage = "Perform a random gauge transform.";
setRandomGauge[] :=
 Do[Block[
   {coords = latticeCoordinates[i],
    vt = randomSUMatrix[]},
  Do[
      setLink[dd, shift[dd, coords, -1],
	      getLink[dd, shift[dd, coords, -1]].vt];
      setLink[dd, coords, ConjugateTranspose[vt].getLink[dd, coords]],
      {dd, nd}]],
    {i, latticeVolume[]}];

setAxialGauge::usage = "Set axial gauge for the current lattice configuration.
That is, links along any Polyakov loop in direction dir are constant and
diagonal.  In the case \"center\" -> True, set the gauge relative to the
nearest element of the center of the group for each axial-direction link.
Option \"abelian\" -> True rotates the field so that it is diagonal.";
Options[setAxialGauge] = {"center" -> False, "abelian" -> True,
                         "debug" -> False};
setAxialGauge[dir_, opts:OptionsPattern[]] :=
 Do[Block[{coords = latticeCoordinates[i], debug = OptionValue["debug"],
       popts = Apply[Sequence, FilterRules[{opts}, Options[getPhases]]]},
   If[coords[[dir]] == 1,
    Block[{v = IdentityMatrix[nc], vt, ct, phases, centerValues},
     (* Initial values *)
     If[False,
	Do[Print[{j, getLink[dir, shift[dir, coords, j - 1]]}],
	   {j, latticeDimensions[[dir]]}]];
     Do[v = v.getLink[dir, shift[dir, coords, j - 1]],
	{j, latticeDimensions[[dir]]}];
     If[debug, Print["Full product:", v]];
     (* Undo "nearest center" that we will add later *)
     If[OptionValue["center"],
        centerValues = Table[
            Last[SUNorm[getLink[dir, shift[dir, coords, j - 1]],
                        "center"->True]],
            {j, latticeDimensions[[dir]]}];
        (* Modify the product to be equivalent to the product
           of the links mod the nearest center element. *)
        v *= Exp[-2 Pi I Total[centerValues]/nc]];
     If[OptionValue["abelian"],
        (* Since we want the diagonal part,
          use getPhases instead of SUPower. *)
        phases = Sort[First[getPhases[v, False]], Greater];
        v = DiagonalMatrix[Exp[I*phases/latticeDimensions[[dir]]]],
        v = SUPower[v, 1/latticeDimensions[[dir]]]];
     If[debug, Print["Root:", v]];
     Do[
      (* Apply gauge choice to each site *)
      ct = shift[dir, coords, j - 1];
      vt = ConjugateTranspose[getLink[dir, shift[dir, ct, -1]]].v;
      If[OptionValue["center"],
         vt *= Exp[I 2 Pi centerValues[[j-1]]/nc]]; 
      Do[setLink[dd, shift[dd, ct, -1],
                 (* In order to decrease the floating point error,
                   just set the axial direction link to the correct value. *)
                 If[dd == dir && !OptionValue["center"],
                    v, getLink[dd, shift[dd, ct, -1]].vt]];
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
Options[setLandauGauge] = {"center" -> False, "directions"->All,
                            "maxAbelianGauge" -> False,
                            "damping" -> 1.7, "debug" -> False};
setLandauGauge::usage = "Set gauge to minimum average link magnitude for the current lattice configuration.  Repeated application should bring the fields into Landau gauge.  The option \"abelian\" -> 0.0 should correspond to maximum abelian gauge.  The option \"directions\" applies the gauge condition only to the specified directions; specifying a single direction would correspond to Axial gauge.";
setLandauGauge[OptionsPattern[]] :=
 Block[{normSum = 0.0, update,
        dirs = If[OptionValue["directions"] === All, nd,
                  OptionValue["directions"]],
        damping = OptionValue["damping"], logLog, ndirs},
   (* Take the log of a link, accumulating statistics
      on the first run through the checkerboard. *)
   logLog = Block[{x = SULog[#1, "center"->OptionValue["center"]], y},
                  If[#2==0, y = Flatten[x];
                            normSum += 2.0*Total[Re[y]^2+Im[y]^2]]; x]&;
   (* Take result from parallel computation and update gaugeField,
      accumulating statistics. *)
   update = (Apply[setLink, Take[#, 3]]; normSum += #[[4]])&;
   ndirs = Sum[1, {dirs}];
   (* Shared variables are super-slow.  Instead, create a table
     of updated links, then update gaugeField and statistics after
     the parallel computation is completed.
     Use checkerboard to avoid conflicts in link updates. *)
   Do[Scan[update, ParallelTable[
       Block[{coords = latticeCoordinates[i], sum, gauge,
              normSum = 0.0},
             If[Mod[Total[coords], 2] == cb,
                (* Construct gauge transform *)
                (* One can use the links themselves
                  and project onto traceless antihermitian
                  matrices, but the performance is a bit worse. *)
                sum = Sum[
                    logLog[getLink[dir, coords], cb]
                    - logLog[getLink[dir, shift[dir, coords, -1]], cb],
                    {dir, dirs}];
                If[OptionValue["maxAbelianGauge"],
                   Do[sum[[i,i]] = 0, {i, nc}]];
                gauge = MatrixExp[-damping*sum/(2*ndirs)];
                (* Create gauge transformed links.
                 Include statistics once. *)
                Table[
                    {{dir, coords, gauge.getLink[dir, coords],
                      If[dir==1, normSum, 0]},
                     {dir, shift[dir, coords, -1],
                      getLink[dir, shift[dir, coords, -1]].
                             ConjugateTranspose[gauge], 0}},
                    {dir, nd}],
                Nothing]],
       {i, latticeVolume[]}, Method->"CoarsestGrained"], {3}], {cb, 0, 1}];
   (* The 2-norm of the initial links,
     averaging over the number of links. *)
   Sqrt[normSum/(latticeVolume[]*nd)]];

setLandauAxialGauge::usage = "Set Axial gauge to minimize norm of all adjoining links.  In practice, this doesn't work very well.";
Options[setMinimumAxialGauge] = {"center" -> False,
                            "damping" -> 1.0, "debug" -> False};
setLaundauAxialGauge[dir1_, OptionsPattern[]] :=
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


nStrategies = 6;
Options[applyGaugeTransforms] = {"abelian" -> False, "maxAbelianGauge"->False,
                                 "loopFunction" -> (Null&)};
applyGaugeTransforms[s_, OptionsPattern[]] :=
 (* New direction each time setAxialGauge[] is called. *)
 Block[{lastDir = nd, dir},
   dir = Function[lastDir = Mod[lastDir, nd] + 1];
   Scan[(Which[
       # == 1, setAxialGauge[dir[], "center" -> False,
                             "abelian"->OptionValue["abelian"]],
       # == 2, setAxialGauge[dir[], "center" -> True,
                             "abelian"->OptionValue["abelian"]],
       # == 3,
       setLandauGauge["center" -> False, "damping" -> 1,
                       "maxAbelianGauge" -> OptionValue["maxAbelianGauge"]],
       # == 4,
       setLandauGauge["center" -> True, "damping" -> 1,
                       "maxAbelianGauge" -> OptionValue["maxAbelianGauge"]],
       # == 5,
       setLandauGauge["center" -> False, "damping" -> 1.5,
                       "maxAbelianGauge" -> OptionValue["maxAbelianGauge"]],
       # == 6,
       setLandauGauge["center" -> True, "damping" -> 1.5,
                       "maxAbelianGauge" -> OptionValue["maxAbelianGauge"]],
       (* In practice, these don't work well *)
       (* # == 7,
        setMinimumAxialGauge[dir[], "center" -> False, "damping" -> 1],
        # == 8,
        setMinimumAxialGauge[dir[], "center" -> True, "damping" -> 1],
        # == 9,
        setMinimumAxialGauge[dir[], "center" -> False, "damping" -> 1.5],
        # == 10,
        setMinimumAxialGauge[dir[], "center" -> True, "damping" -> 1.5]*)
       True,
       Abort[]
         ];OptionValue["loopFunction"][])&, s];
   latticeNorm["center" -> True]];

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
             Less[#1[[1, 2]], #2[[1, 2]]]&]];


dualShift3 = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
makeFaradayLattice::usage = 
  "Construct a dual lattice containing a^2*g times the Faraday tensor.\
This construction only makes sense in the context of Landau gauge.";
makeFaradayLattice[] := 
 Block[{u1, u2, u1s2, u2s1, faces = {{2, 3}, {3, 1}, {1, 2}}, 
    dualLattice = Array[Null &, {nd, latticeVolume[]}], dir1, dir2, 
    coords, z}, 
  Do[{dir1, dir2} = faces[[dir0]]; 
      coords = latticeCoordinates[i];
    u1 = getLink[dir1, coords];
    u2s1 = getLink[dir2, shift[dir1, coords]];
    u1s2 = getLink[dir1, shift[dir2, coords]];
    u2 = getLink[dir2, coords];
    coords = wrapIt[coords + dualShift3[[dir0]]];
    z = u1.u2s1 + ConjugateTranspose[u2.u1s2] +
         u2s1.ConjugateTranspose[u1s2] + ConjugateTranspose[u2].u1;
    dualLattice[[dir0, linearSiteIndex[coords]]] =
    -I (z - ConjugateTranspose[z])/4,
     {dir0, nd}, {i, latticeVolume[]}]; dualLattice]/;nd==3;

makeFaradayLatticeA::usage = "Alternative form for makeFaradayLattice[]."; 
makeFaradayLatticeA[testTerm_:1.0] := 
 Block[{lgf = (-I)*ParallelMap[SULog, gaugeField, {2}], agA1, agA2, agA1s2, 
         agA2s1, agA11, agA22, faces = {{2, 3}, {3, 1}, {1, 2}}, 
    dualLattice = Array[Null &, {nd, latticeVolume[]}], dir1, dir2, 
    coords},
  Do[{dir1, dir2} = faces[[dir0]];
      coords = latticeCoordinates[i];
    agA1 = getLink[lgf][dir1, coords];
    agA2s1 = getLink[lgf][dir2, shift[dir1, coords]];
    agA1s2 = getLink[lgf][dir1, shift[dir2, coords]];
    agA2 = getLink[lgf][dir2, coords];
    agA11 = (agA1 + agA1s2)/2; agA22 = (agA2 + agA2s1)/2;
    coords = wrapIt[coords + dualShift3[[dir0]]];
    dualLattice[[dir0, linearSiteIndex[coords]]] =
    (agA2s1 - agA2) - (agA1s2 - agA1) +
    testTerm*I*(agA11.agA22 - agA22.agA11),
     {dir0, nd}, {i, latticeVolume[]}]; dualLattice]/;nd==3;
