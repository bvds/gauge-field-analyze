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
latticeCoordinates[i_] := latticeCoordinates[i, latticeDimensions];
latticeCoordinates[ii_, dimensions_] :=
 Block[{i = ii - 1, j},
  Map[({i, j} = QuotientRemainder[i, #]; j + 1)&, dimensions]];

latticeIndex::usage = "Inverse of latticeCoordinates.";
latticeIndex[coords_List] := latticeIndex[coords, latticeDimensions];
latticeIndex[coords_List, dimensions_] :=
 Block[{i = 1, delta = 1},
  Do[i += (coords[[k]] - 1)*delta;
   delta *= dimensions[[k]], {k, Length[dimensions]}]; i];

makeCoordList[] := Block[
    {result = Array[Null&, latticeVolume[]]},
    Do[Block[{coord = latticeCoordinates[k]},
             result[[linearSiteIndex[coord]]] = coord],
       {k, latticeVolume[]}]; result];

shift::usage = "Wrap-around, assuming periodic boundary conditions.";
shift[dir_, coordIn_List, size_:1] :=
  Block[{cc = coordIn},
	cc[[dir]] = 1 + Mod[cc[[dir]] + size - 1, latticeDimensions[[dir]]];
	cc];

wrapIt[coords_] := wrapIt[coords, latticeDimensions];
wrapIt[coords_List, dims_] :=
    MapThread[(1 + Mod[#1 - 1, #2])&, {coords, dims}];

getLatticeSlice::usage = "Construct a slice of the lattice with the given direction and coordinate value.  Returns {latticeDimension, gaugeField} using latticeIndex[coord] ordering.";
getLatticeSlice[dir0_, t_] :=
    Block[
        {dims = Drop[latticeDimensions, {dir0}]},
        {dims, Table[
            Block[{coord = latticeCoordinates[k, dims]},
                  getLink[dir1, Insert[coord, t, dir0]]],
            {dir1, Drop[Range[nd], {dir0}]}, {k, latticeVolume[dims]}]}
    ];


(* Handle lattice where some number of links-sites have
  been frozen. *)

lattice::invalidBC = "Unknown value for latticeBC"
boundaryLinkQ[dir_, coords_] :=
    Which[latticeBC == "PERIODIC_GAUGEBC",
          False,
          latticeBC == "TRANS_GAUGEBC",
          Block[{f = ((dir != # && coords[[#]] == 1) ||
                      coords[[#]] == latticeDimensions[[#]])&},
                f[transverseBCDirections[[1]]] ||
                f[transverseBCDirections[[2]]]],
          True,
          Message[lattice::ivalidBC];
          $Failed];
boundarySiteQ[coords_] :=
    Which[latticeBC == "PERIODIC_GAUGEBC",
          False,
          latticeBC == "TRANS_GAUGEBC",
          Block[{f = (coords[[#]] == 1 || coords[[#]] ==
                      latticeDimensions[[#]])&},
                f[transverseBCDirections[[1]]] ||
                f[transverseBCDirections[[2]]]],
          True,
          Message[lattice::ivalidBC];
          $Failed];
reducedSiteCount[] :=
    latticeVolume[]*Which[
        latticeBC == "PERIODIC_GAUGEBC",
        1,
        latticeBC == "TRANS_GAUGEBC",
        Block[{l1, l2},
              {l1, l2} = Map[latticeDimensions[[#]]&,
                             transverseBCDirections];
              (l1-2)*(l2-2)/(l1 l2)],
        True,
        Message[lattice::invalidBC];
        $Failed];
reducedLinkCount[] :=
    latticeVolume[]*Which[
        latticeBC == "PERIODIC_GAUGEBC",
        nd,
        latticeBC == "TRANS_GAUGEBC",
        Block[{l1, l2},
              {l1, l2} = Map[latticeDimensions[[#]]&,
                             transverseBCDirections];
              (nd*(l1-2)*(l2-2) + l1 + l2 - 4)/(l1 l2)],
        True,
        Message[lattice::invalidBC];
        $Failed];
reducedSiteIndex::usage = "Version of latticeIndex where only sites that have not been masked are indexed.";
reducedSiteIndex[coords_] :=
    Which[
        latticeBC == "PERIODIC_GAUGEBC",
        latticeIndex[coords],
        latticeBC == "TRANS_GAUGEBC",
        Block[{cords = coords, dims = latticeDimensions},
              Do[
                  cords[[transverseBCDirections[[i]]]] -= 1;
                  dims[[transverseBCDirections[[i]]]] -= 2, {i, 2}];
              latticeIndex[cords, dims]],
        True,
        Message[lattice::invalidBC];
        $Failed];
(* Helper function to construct a lookup table since
  a direct calculation is too difficult. *)
reducedLinkIndex::count = "Wrong number of links. cache length `1`, count=`2`, reducedLinkCount = `3`";
reducedLinkCache[dims_List, transDirs_List] :=
    reducedLinkCache[dims, transDirs] =
    Block[{y = Association[], count = 1},
          Do[
              Block[{coords = latticeCoordinates[i]},
                    If[!boundaryLinkQ[dir, coords],
                       y[{dir, coords}] = count++]],
              {i, latticeVolume[]}, {dir, nd}];
          If[Length[y] != reducedLinkCount[],
             (* Print[{dims, transDirs, latticeDimensions}]; *)
             Message[reducedLinkIndex::count,
                     Length[y], count, reducedLinkCount[]]];
          y];
reducedLinkIndex[dir_, coords_List] :=
    Which[
        latticeBC == "PERIODIC_GAUGEBC",
        (latticeIndex[coords] -1) * nd + dir,
        latticeBC == "TRANS_GAUGEBC",
        reducedLinkCache[latticeDimensions, transverseBCDirections][
            {dir, coords}],
        True,
        Message[lattice::invalidBC];
        $Failed];


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
  "All links identity for a given nd, nc, latticeDimensions.
Can add small random pertubation, optionally diagonal,  or choose \
a random element of the center.";
Options[makeTrivialLattice] = {randomGauge -> False,
  randomCenter -> False, "diagonal" -> False,
  randomPerturbation -> 0};
makeTrivialLattice[OptionsPattern[]] := (
    gaugeField =
    Table[If[
        OptionValue[randomPerturbation] > 0,
        MatrixExp[I If[OptionValue["diagonal"],
                       RandomReal[{-1, 1}*OptionValue[randomPerturbation],
	                          nc - 1].Take[SUGenerators[], -(nc - 1)],
                       RandomReal[{-1, 1}*OptionValue[randomPerturbation],
	                          nc^2 - 1].SUGenerators[]]],
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


(* See files "maxwell-lattice.nb" and "coulomb.nb"
  Second argument is length of wire, third is separation. *)
Clear[wilsonCoulomb];
wilsonCoulomb[False, l_, r_, eps_:1] := wilsonCoulomb[False, l, r, eps] =
    Block[{eta = Sqrt[l^2 + r^2]},
          2*(l*Log[l/eps] - l - r + eta - l*Log[(eta + l)/r])];
wilsonCoulomb[True, l_?NumericQ, r_?NumericQ] := wilsonCoulomb[True, l, r] =
Sum[pointPotential3[l1 - l2, 0, 0] -
    pointPotential3[l1 - l2, r, 0], {l1, l}, {l2, l}];


(* For these constants to be used in plots and tables, they
  are defined as global variables.
  Alternatively, one could wrap these and the function in Module[]
  to define a unique instance. *)
Clear[chiSquared, a2sigma, cNorm, logNorm, c0, c1, c2, c3, c4, c5, c6, cm1,
      coff1, coff3, coff4,
      cCoulomb, cPerimeter,
      eigenNorm, sConstant, lConstant, rConstant, offConstant];
Format[chiSquared] = Superscript["\[Chi]", 2];
Format[a2sigma] = Row[{Style["a", Italic]^2, "\[Sigma]"}];
Format[cNorm] = Subscript["c", "z"];
Format[cNorm[i_]] = Row[cNorm, "(", i, ")"];
Format[logNorm] = Log[cNorm];
Format[logNorm[i_]] = Log[cNorm[i]];
Format[cCoulomb] = Subscript["c", "q"];
Format[cPerimeter] = Subscript["c", "p"];
Do[Format[Symbol["c"<>ToString[i]]] = Subscript["c", i], {i, 0, 6}];
Do[Format[Symbol["cm"<>ToString[i]]] = Subscript["c", -i], {i, 1}];
Format[eigenNorm[i_]] := Subscript["C", i];
Format[sConstant[i_]] := Row[{Style["a", Italic]^2,
                              Subscript["\[Sigma]", i]}];
Format[lConstant[i_]] := Subscript["d", i];
Format[rConstant[i_]] := Subscript["e", i];
Format[rConstant[i_, j_]] := Subscript["e", Row[{i, j}]];
Format[offConstant[i_]] := Subscript["f", i];
Format[offConstant[i_, j_]] := Subscript["f", Row[{i, j}]];
rescaleCovarianceMatrix::usage = "This was mostly experimental; see comments in \"gauge.nb\"";
rescaleCovarianceMatrix[cov_?MatrixQ, err_?VectorQ, label_: None] :=
 Block[{rescale = err/Sqrt[Diagonal[cov]]}, 
  Print["Rescale:  ", {label, Exp[ Mean[Log[rescale]]], 
    Exp[StandardDeviation[Log[rescale]]]}]; 
  Map[(#*rescale) &, cov*rescale]]

casimir::usage = "For the Nambu-Goto string picture.  Since we average
        over all pairs of a given separation, this will only couple
        to the even parity states.";
casimir[state_] :=
    Block[{degeneracy = Function[i, {0, 1, 2, 2}[[i + 1]]]},
          8 Pi (degeneracy[state] - (nd-2)/24)];


stringModelBare::usage = "Fit to an area law plus power law corrections, including the universal string corrections.  See Andreas Athenodorou, Barak Bringoltz, Michael Teper https://arxiv.org/abs/0709.0693; Ofer Aharony & Zohar Komargodski arXiv:1302.6257v2 [hep-th] 12 Mar 2013; Teper review article http://arxiv.org/abs/0912.3339  For option \"stringTension\"->Automatic (default), fit string tension, else use the value given.  Likewise for options \"coulomb\" and \"perimeter.\"";
Options[stringModelBare] = {printResult -> False, "lowerCutoff" -> 1/2,
                        "upperCutoff" -> Infinity, "order" -> 0,
                        "stringCorrections" -> False,
                        "eigenstate" -> 0, "NambuGoto" -> False,
                        "offAxis" -> 0, "linearLCorrections" -> True,
                        "pointPotential" -> True,
                        "rescaleCovarianceMatrix" -> False,
                        "stringTension" -> Automatic,
                        "coulomb" -> 0, "perimeter" -> Automatic,
                        "covarianceMatrix" -> None};
stringModel::lowerCutoff = "Invalid combination of \"pointPotential\" and \"lowerCutoff,\" exiting";
stringModelBare[tallyData_, OptionsPattern[]] :=
 Block[(* Protect against any global definitions of model parameters
         and allow for local constant value. *)
     {ff, f0, x, y, ll, logNorm, a2sigma, cCoulomb, cPerimeter,
      c0, c1, c2, c3, c4, c5, 
      offConstant, rConstant, lConstant,
      coff1, coff3, coff4,
      potentialForm = (OptionValue["eigenstate"] == 0),
      zeroQ = Function[x, NumericQ[x] && x==0],
      nll = Length[Union[Map[Last, Keys[tallyData]]]],
      llValues = Union[Map[Last, Keys[tallyData]]],
      na = Length[First[Keys[tallyData]]] - 1,
      data, cov,
      sin2theta2 = Function[{x, y}, (2 x y/(x^2 + y^2))^2],
      filter = Map[First, Select[
          MapIndexed[{#2[[1]], #1}&, Keys[tallyData]],
          (Norm[#[[2,1]]] > OptionValue["lowerCutoff"] &&
           Norm[#[[2,1]]] < OptionValue["upperCutoff"])&]]},
  data = Normal[tallyData][[filter]];
  cov = If[OptionValue["covarianceMatrix"] === None,
           Map[(#[[2, 2]]^2)&, data],
           Block[{cc = OptionValue["covarianceMatrix"][[filter, filter]]},
                 If[OptionValue["rescaleCovarianceMatrix"],
                    rescaleCovarianceMatrix[cc, Map[#[[2,2]]&, data]],
                    cc]]];
  If[Not[OptionValue["pointPotential"] || OptionValue["lowerCutoff"] > 0],
     Message[stringModel::lowerCutoff];
     Return[$Failed]];
  If[OptionValue["coulomb"] =!= Automatic,
     cCoulomb = OptionValue["coulomb"]];
  Which[OptionValue["perimeter"] =!= Automatic,
     cPerimeter = OptionValue["perimeter"],
     (* With only one longitudinal value, this can
       be absorbed into logNorm. *)
     nll<2, cPerimeter = 0];
  If[OptionValue["stringTension"] =!= Automatic,
     a2sigma = OptionValue["stringTension"]];
  (* With only one longitudinal value, this is indistinguishable
    from the string tension.
    If nll==1, assume cPerimeter would use up L-dependence. *)
  Which[OptionValue["order"]<2 || nll<3, c0 = 0,
        (* Universal closed string leading correction *)
        OptionValue["stringCorrections"], c0 = - N[Pi (nd - 2)/8]];
  Which[OptionValue["order"]<2, c1 = 0,
        (* Universal open string leading correction *)
        OptionValue["stringCorrections"], c1 = - N[Pi (nd - 2)/24]];
  If[OptionValue["order"]<3, c3 = 0];
  If[OptionValue["order"]<3 || !OptionValue["linearLCorrections"], c4 = 0];
  (* With three longitudinal values, this can be
    absorbed into logNorm, cPerimeter, and c2. *)
  If[OptionValue["order"]<3 || nll<4, c5 = 0];
  If[OptionValue["offAxis"]<2, coff1 = 0];
  If[OptionValue["offAxis"]<3, coff3 = 0];
  If[OptionValue["offAxis"]<3 || !OptionValue["linearLCorrections"], coff4 = 0];
  Do[If[OptionValue["order"]<i, rConstant[j, i] = 0];
     If[OptionValue["offAxis"]<i, offConstant[j, i] = 0],
     {j, 0, OptionValue["eigenstate"] - 1}, {i, 1, 2}];
  (* Perform a simple fit to get decent starting values *)
  f0 = covariantFit2[
      cov,
      Map[{Norm[#[[1, 1]]], #[[1,-1]], #[[2, 1]]}&, data],
      Exp[-r*ll*a2sigma + logNorm],
      {{logNorm, 0.05},
       If[OptionValue["stringTension"] === Automatic,
          {a2sigma, 0.01888}, Nothing]},
      {r, ll},
      Method -> "LevenbergMarquardt"];
  ff = covariantFit2[
      cov,
      Map[Append[Flatten[#[[1]]], #[[2, 1]]]&, data],
      (* The constrained version of NonLinearModelFit[] does not  
        use Levenberg-Marquardt and does not converge as well.
        Instead, add boundary terms to the fit function.

        Tried Exp[-const*r[i]] factor in Coulomb term for 12x18x18,
        but preferred sign was const<0. *)
      Sum[Block[
          {r = Sqrt[x[i]^2 + y[i]^2],
           offAxis = sin2theta2[x[i], y[i]],
           pointPot = OptionValue["pointPotential"]}, 
          Which[
              potentialForm,
              (* area law potential shape with Coulomb correction
                and power law  corrections *)
              Exp[-r*ll*a2sigma + logNorm -
                     cCoulomb*ll*If[pointPot,
                         pointPotential2[x[i], y[i]], Log[r]] -
                     2*cPerimeter*ll -
                     (* second order corrections *)
                     (* Luscher term, resummed *)
                     c0*r/ll -
                     (* leading r correction to string tension *)
                     (c1 + coff1*offAxis)*ll/r -
                     (* third order corrections *)
                     (* leading r correction to coulomb potential *)
                     (c3 + coff3*offAxis)*ll/r^2 -
                     (* linear correction, maybe not allowed? *)
                     (c4 + coff4*offAxis)/r -
                     (* second Luscher term, resummed *)
                     c5/ll],
              (* Match Nambu-Goto spectrum.  Ignore the small
                L correction, even for excited states, since we
                have other, larger errors. *)
              OptionValue["NambuGoto"],
              Sum[eigenNorm[j]*Exp[
                  -r*Sqrt[Max[
                      (ll*a2sigma)^2 + a2sigma*casimir[j], 10^-10]] -
                  (rConstant[j, 1] + offConstant[j, 1]*offAxis)*ll/r -
                  (rConstant[j, 2] + offConstant[j, 2]*offAxis)*ll/r^2],
                  {j, 0, OptionValue["eigenstate"] - 1}],
              True,
              (* Fit spectrum, but don't make any assumption
                about agreement with Nambu-Goto.  However,
                assume the leading correction is O(1/L). *)
              Sum[Exp[eigenNorm[j]
                  -r*ll*If[j==0, a2sigma, sConstant[j]] -
                  lConstant[j]*r/ll -
                  (rConstant[j, 1] - offConstant[j, 1]*offAxis)*ll/r -
                  (rConstant[j, 2] - offConstant[j, 2]*offAxis)*ll/r^2],
                  {j, 0, OptionValue["eigenstate"] - 1}]
          ]], {i, na}],
      Join[
          {If[!NumericQ[a2sigma],
              {a2sigma, a2sigma/.f0["BestFitParameters"]},  Nothing]},
          If[potentialForm,
             {{logNorm, logNorm/.f0["BestFitParameters"]},
              If[!NumericQ[cCoulomb], {cCoulomb, 0.02}, Nothing],
              If[!NumericQ[cPerimeter], {cPerimeter, 0.}, Nothing],
              If[!NumericQ[c0], {c0, 0.0}, Nothing],
              If[!NumericQ[c1], {c1, 0.0}, Nothing],
              If[!NumericQ[c3], {c3, 0.0}, Nothing],
              If[!NumericQ[c4], {c4, 0.0}, Nothing],
              If[!NumericQ[c5], {c5, 0.0}, Nothing],
              If[!NumericQ[coff1], {coff1, 0.0}, Nothing],
              If[!NumericQ[coff3], {coff3, 0.0}, Nothing],
              If[!NumericQ[coff4], {coff4, 0.0}, Nothing]},
             Join[
                 Table[{eigenNorm[j], logNorm/(j+1.0)/.f0["BestFitParameters"]},
                       {j, 0, OptionValue["eigenstate"] - 1}],
                 If[!OptionValue["NambuGoto"],
                    Table[{sConstant[j], (j+1)*
                          If[NumericQ[a2sigma],
                             a2sigma, a2sigma/.f0["BestFitParameters"]]},
                          {j, 1, OptionValue["eigenstate"] - 1}],
                    {}],
                 Table[If[!NumericQ[lConstant[j]] && !OptionValue["NambuGoto"],
                          {lConstant[j], 0.0}, Nothing],
                       {j, 0, OptionValue["eigenstate"] - 1}],
                 Flatten[Table[
                     If[!NumericQ[rConstant[j, i]],
                        {rConstant[j, i], 0.0}, Nothing],
                     {j, 0, OptionValue["eigenstate"] -1}, {i, 2}], 1],
                 Flatten[Table[
                     If[!NumericQ[offConstant[j, i]],
                        {offConstant[j, i], 0.0}, Nothing],
                     {j, 0, OptionValue["eigenstate"] -1}, {i, 2}], 1]
             ]]],
      Append[Flatten[Table[{x[i], y[i]}, {i, na}]], ll],
      Method -> "LevenbergMarquardt"];
  ff["filter"] = filter;
  If[OptionValue[printResult],
     Print["Correlation matrix: ",ff["CorrelationMatrix"]];
     (* Print["Covariance matrix: ", ff["CovarianceMatrix"]]; *)
     Print[ff["ParameterTable"]];
     Print["chi^2: ", ff["chiSquared"],
           " for ",
           Length[data]-Length[ff["BestFitParameters"]], " d.o.f."]];
  ff];


stringModelBlocked::usage = "For blocked operators, fit to an area law plus power law corrections, including the universal string corrections.  See Andreas Athenodorou, Barak Bringoltz, Michael Teper https://arxiv.org/abs/0709.0693; Ofer Aharony & Zohar Komargodski arXiv:1302.6257v2 [hep-th] 12 Mar 2013; Teper review article http://arxiv.org/abs/0912.3339.  Use a cosh for periodicity, https://arxiv.org/abs/hep-lat/0107007.  For option \"stringTension\"->Automatic (default), fit string tension, else use the value given.  Likewise for options \"coulomb\" and \"perimeter.\"";
Options[stringModelBlocked] = {printResult -> False, "lowerCutoff" -> 1/2,
                        "upperCutoff" -> Infinity, "order" -> 0,
                        "stringCorrections" -> False,
                        "eigenstate" -> 0, "NambuGoto" -> False,
                        "linearLCorrections" -> True,
                        "pointPotential" -> True,
                        "rescaleCovarianceMatrix" -> False,
                        "stringTension" -> Automatic,
                        (* Since we are blocking differently for
                          each L (lattice size in the direction
                          of the Polyakov loop), the normalization will
                          normally be different for each L, making the
                          perimeter term indistinguishable.*)
                        "longitudinalNorms" -> True,
                        "coulomb" -> 0, "perimeter" -> 0,
                        "covarianceMatrix" -> None};
stringModelBlocked[tallyData_, OptionsPattern[]] :=
 Block[(* Protect against any global definitions of model parameters
         and allow for local constant value. *)
     {ff, f0, ll, logNorm, a2sigma, cCoulomb, cPerimeter,
      c0, c1, c2, c3, c4, c5, 
      rConstant, lConstant,
      potentialForm = (OptionValue["eigenstate"] == 0),
      nll = Length[Union[Map[Last, Keys[tallyData]]]],
      llValues = Union[Map[Last, Keys[tallyData]]],
      na = Length[First[Keys[tallyData]]] - 1,
      data, cov,
      filter = Map[First, Select[
          MapIndexed[{#2[[1]], #1}&, Keys[tallyData]],
          (#[[2,1]] > OptionValue["lowerCutoff"] &&
           #[[2,1]] < OptionValue["upperCutoff"])&]]},
  data = Normal[tallyData][[filter]];
  cov = If[OptionValue["covarianceMatrix"] === None,
           Map[(#[[2, 2]]^2)&, data],
           Block[{cc = OptionValue["covarianceMatrix"][[filter, filter]]},
                 If[OptionValue["rescaleCovarianceMatrix"],
                    rescaleCovarianceMatrix[cc, Map[#[[2,2]]&, data]],
                    cc]]];
  If[Not[OptionValue["pointPotential"] || OptionValue["lowerCutoff"] > 0],
     Message[stringModel::lowerCutoff];
     Return[$Failed]];
  If[OptionValue["coulomb"] =!= Automatic,
     cCoulomb = OptionValue["coulomb"]];
  Which[OptionValue["perimeter"] =!= Automatic,
     cPerimeter = OptionValue["perimeter"],
     (* With only one longitudinal value, this can
       be absorbed into logNorm. *)
     nll<2, cPerimeter = 0];
  If[OptionValue["stringTension"] =!= Automatic,
     a2sigma = OptionValue["stringTension"]];
  (* With only one longitudinal value, this is indistinguishable
    from the string tension.
    If nll==1, assume cPerimeter would use up L-dependence. *)
  Which[OptionValue["order"]<2 || nll<3, c0 = 0,
        (* Universal closed string leading correction *)
        OptionValue["stringCorrections"], c0 = - N[Pi (nd - 2)/8]];
  Which[OptionValue["order"]<2, c1 = 0,
        (* Universal open string leading correction *)
        OptionValue["stringCorrections"], c1 = - N[Pi (nd - 2)/24]];
  If[OptionValue["order"]<3, c3 = 0];
  If[OptionValue["order"]<3 || !OptionValue["linearLCorrections"], c4 = 0];
  (* With three longitudinal values, this can be
    absorbed into logNorm, cPerimeter, and c2. *)
  If[OptionValue["order"]<3 || nll<4, c5 = 0];
  If[OptionValue["offAxis"]<3 || !OptionValue["linearLCorrections"], coff4 = 0];
  Do[If[OptionValue["order"]<i, rConstant[j, i] = 0],
     {j, 0, OptionValue["eigenstate"] - 1}, {i, 1, 2}];
  (* Perform a simple fit to get decent starting values *)
  f0 = covariantFit2[
      cov,
      Map[Append[#[[1]], #[[2, 1]]]&, data],
      Exp[-(r-rr/2)*ll*a2sigma + logNorm],
      {{logNorm, 0.1},
       If[OptionValue["stringTension"] === Automatic,
          {a2sigma, 0.01888}, Nothing]},
      {r, rr, ll},
      Method -> "LevenbergMarquardt"];
  ff = covariantFit2[
      cov,
      Map[Append[#[[1]], #[[2, 1]]]&, data],
      (* The constrained version of NonLinearModelFit[] does not  
        use Levenberg-Marquardt and does not converge as well.
        Instead, add boundary terms to the fit function.

        Tried Exp[-const*r[i]] factor in Coulomb term for 12x18x18,
        but preferred sign was const<0. *)
      Sum[Block[
          {pointPot = OptionValue["pointPotential"]}, 
          Which[
              potentialForm,
              (* area law potential shape with Coulomb correction
                and power law  corrections *)
              Exp[-r*ll*a2sigma + If[OptionValue["longitudinalNorms"],
                                     logNorm[ll], logNorm] -
                     cCoulomb*ll*If[pointPot,
                         pointPotential2[x[i], y[i]], Log[r]] -
                     2*cPerimeter*ll -
                     (* second order corrections *)
                     (* Luscher term, resummed *)
                     c0*r/ll -
                     (* leading r correction to string tension *)
                     c1*ll/r -
                     (* third order corrections *)
                     (* leading r correction to coulomb potential *)
                     c3*ll/r^2 -
                     (* linear correction, maybe not allowed? *)
                     c4/r -
                     (* second Luscher term, resummed *)
                     c5/ll],
              (* Match Nambu-Goto spectrum.  Ignore the small
                L correction, even for excited states, since we
                have other, larger errors. *)
              OptionValue["NambuGoto"],
              Sum[eigenNorm[j]*Exp[
                  -r*Sqrt[Max[
                      (ll*a2sigma)^2 + a2sigma*casimir[j], 10^-10]] -
                  rConstant[j, 1]*ll/r -
                  rConstant[j, 2]*ll/r^2],
                  {j, 0, OptionValue["eigenstate"] - 1}],
              True,
              (* Fit spectrum, but don't make any assumption
                about agreement with Nambu-Goto.  However,
                assume the leading correction is O(1/L). *)
              Sum[Exp[eigenNorm[j]
                  -r*ll*If[j==0, a2sigma, sConstant[j]] -
                  lConstant[j]*r/ll -
                  rConstant[j, 1]*ll/r -
                  rConstant[j, 2]*ll/r^2],
                  {j, 0, OptionValue["eigenstate"] - 1}]
          ]], {r, {r0, rr - r0}}],
      Join[
          {If[!NumericQ[a2sigma],
              {a2sigma, a2sigma/.f0["BestFitParameters"]},  Nothing]},
          If[potentialForm,
             Join[
                 If["longitudinalNorms",
                    Table[{logNorm[i], logNorm/.f0["BestFitParameters"]},
                          {i, nl}],
                    {{logNorm, logNorm/.f0["BestFitParameters"]}}],
                 {If[!NumericQ[cCoulomb], {cCoulomb, 0.02}, Nothing],
                  If[!NumericQ[cPerimeter], {cPerimeter, 0.}, Nothing],
                  If[!NumericQ[c0], {c0, 0.0}, Nothing],
                  If[!NumericQ[c1], {c1, 0.0}, Nothing],
                  If[!NumericQ[c3], {c3, 0.0}, Nothing],
                  If[!NumericQ[c4], {c4, 0.0}, Nothing],
                  If[!NumericQ[c5], {c5, 0.0}, Nothing]}],
             Join[
                 Table[{eigenNorm[j], Exp[logNorm]/(j+1.0)/.f0["BestFitParameters"]},
                       {j, 0, OptionValue["eigenstate"] - 1}],
                 If[!OptionValue["NambuGoto"],
                    Table[{sConstant[j], (j+1)*
                          If[NumericQ[a2sigma],
                             a2sigma, a2sigma/.f0["BestFitParameters"]]},
                          {j, 1, OptionValue["eigenstate"] - 1}],
                    {}],
                 Table[If[!NumericQ[lConstant[j]] && !OptionValue["NambuGoto"],
                          {lConstant[j], 0.0}, Nothing],
                       {j, 0, OptionValue["eigenstate"] - 1}],
                 Flatten[Table[
                     If[!NumericQ[rConstant[j, i]],
                        {rConstant[j, i], 0.0}, Nothing],
                     {j, 0, OptionValue["eigenstate"] -1}, {i, 2}], 1]]]],
      {r0, rr, ll},
      Method -> "LevenbergMarquardt"];
  ff["filter"] = filter;
  If[OptionValue[printResult],
     Print["Correlation matrix: ",ff["CorrelationMatrix"]];
     (* Print["Covariance matrix: ", ff["CovarianceMatrix"]]; *)
     Print[ff["ParameterTable"]];
     Print["chi^2: ", ff["chiSquared"],
           " for ",
           Length[data]-Length[ff["BestFitParameters"]], " d.o.f."]];
  ff];

wilsonModel::usage = "Fit to an exponential, including Coulomb force contributions and an optional excited state.  The Coulomb force contributions consist of a complicated normal term plus a perimeter term.  For option \"stringTension\"->Automatic (default), fit string tension, else use value given.  Likewise for \"coulomb\" and \"perimeter.\"";
Options[wilsonModel] = {printResult -> False, "order" -> 0,
                        (* lowerCutoff limits the length of the shortest side;
                          upperCutoff limits the lattice size minus the width
                          of the Wilson loop. *)
                        "lowerCutoff" -> 0, "upperCutoff" -> 0,
                        "stringTension" -> Automatic,
                        "stringCorrections" -> False,
                        "pointPotential" -> True,
                        "coulomb" -> 0, "perimeter" -> Automatic,
                        "covarianceMatrix" -> None};
wilsonModel[data0_, OptionsPattern[]] :=
 Block[(* Protect against any global definitions of model parameters
         and allow for local constant value. *)
     {ff, w1, w2, l1, l2,
      a2sigma, cCoulomb, cPerimeter, cm1, c1, c2, c3, c4, c5, c6,
      data,
     filter = Map[First, Select[
          MapIndexed[{#2[[1]], #1}&, Keys[data0]],
          (#[[2, 1]] > OptionValue["lowerCutoff"] &&
           #[[2, 2]] > OptionValue["lowerCutoff"] &&
           #[[2, 1]] < #[[2, 3]] - OptionValue["upperCutoff"] &&
           #[[2, 2]] < #[[2, 4]] - OptionValue["upperCutoff"])&]]},
  data = Normal[data0][[filter]];
  If[OptionValue["stringTension"] =!= Automatic,
     a2sigma = OptionValue["stringTension"]];
  If[OptionValue["coulomb"] =!= Automatic,
     cCoulomb = OptionValue["coulomb"],
     If[OptionValue["order"]<1, cCoulomb = 0]];
  If[OptionValue["perimeter"] =!= Automatic,
     cPerimeter = OptionValue["perimeter"],
     If[OptionValue["order"]<1, cPerimeter = 0]];
  (* c1 = 0; c5=0; c3 = 0; c6 = 0; *)
  If[OptionValue["order"]<1, cm1 = 0];
  Which[OptionValue["order"]<1, c1 = 0,
        (* Universal open string leading correction *)
        OptionValue["stringCorrections"], c1 = - N[Pi (nd - 2)/24]];
  If[OptionValue["order"]<2, c2 = 0];
  Which[OptionValue["order"]<2, c3 = 0,
        (* This term is universally zero for all strings *)
        OptionValue["stringCorrections"], c3 = 0];
  If[OptionValue["order"]<3, c4 = 0];
  If[OptionValue["order"]<3, c5 = 0];
  Which[OptionValue["order"]<3, c6 = 0,
     (* This term is also universal and small, 0.00857 (D-2) *)
     OptionValue["stringCorrections"], c6 = -N[Pi^2 (nd - 2)^2/1152]];
  ff = covariantFit2[
    If[OptionValue["covarianceMatrix"] === None,
       Map[(#[[2, 2]]^2)&, data],
       OptionValue["covarianceMatrix"][[filter, filter]]],
    Map[Append[#[[1]], #[[2, 1]]] &, data], 
    Exp[logNorm-w1*w2*a2sigma -
           cCoulomb*(wilsonCoulomb[OptionValue["pointPotential"], w1, w2] +
                     wilsonCoulomb[OptionValue["pointPotential"], w2, w1]) -
           cPerimeter*2*(w1 + w2) -
           c1*(w1/w2 + w2/w1) -
           c2*(1/w1 + 1/w2) -
           c3*(w1/w2^2 + w2/w1^2) -
           c4/(w1*w2) -
           c5*(1/w1^2 + 1/w2^2) -
           c6*(w1/w2^3 + w2/w1^3)] +
    (* Contribution from the exterior of the loop,
      using lattice periodicity.

      Dividing into two terms, one for each direction,
      gives an extremely poor fit. *)
    Exp[logNorm-(l1*l2 - w1*w2)*a2sigma -
           cCoulomb*(wilsonCoulomb[OptionValue["pointPotential"],
                                   w1, l2 - w2] +
                     wilsonCoulomb[OptionValue["pointPotential"],
                                   w2, l1 - w1]) -
           cPerimeter*2*(w1 + w2) - cm1*(l1 + l2 - w1 - w2)],
    {If[!NumericQ[a2sigma], {a2sigma, 0.016}, Nothing],
     {logNorm, 1.0},
     If[!NumericQ[cCoulomb], {cCoulomb, 0.01}, Nothing],
     If[!NumericQ[cPerimeter], {cPerimeter, 0.01}, Nothing],
     If[!NumericQ[cm1], {cm1, 0.0}, Nothing],
     If[!NumericQ[c1], {c1, 0.0}, Nothing],
     If[!NumericQ[c2], {c2, 0.0}, Nothing],
     If[!NumericQ[c3], {c3, 0.0}, Nothing],
     If[!NumericQ[c4], {c4, 0.0}, Nothing],
     If[!NumericQ[c5], {c5, 0.0}, Nothing],
     If[!NumericQ[c6], {c6, 0.0}, Nothing]},
    {w1, w2, l1, l2},
    Method -> "LevenbergMarquardt"];
  If[OptionValue[printResult],
     Print["Correlation matrix: ",ff["CorrelationMatrix"]];
     (* Print["Covariance matrix: ",ff["CovarianceMatrix"]]; *)
     Print[ff["ParameterTable"]];
     Print[ColumnForm[Map[
                  Row[If[#[[1,1]] === logNorm,
                         {cNorm, " = ",
                          Exp[valueError[#[[1, 2]], #[[2]]]]},
                         {#[[1, 1]], " = ", 
                          valueError[#[[1, 2]], #[[2]]]}]]&, 
            Transpose[
                ff[{"BestFitParameters", 
                    "ParameterErrors"}]]]]];
     Print["chi^2: ", ff["chiSquared"],
           " for ",
           Length[data]-Length[ff["BestFitParameters"]], " d.o.f."]];
  ff];


setGlobalGauge::usage = "Perform global color rotation. u is a unitary matrix.";
setGlobalGauge[u_] :=
    (gaugeField = Map[u.#.ConjugateTranspose[u]&, gaugeField, {2}]);

setSiteGauge::usage == "Apply a gauge transform to a single lattice site.";
setSiteGauge[coords_, u_] := 
    Do[setLink[dir, coords, u.getLink[dir, coords]];
       setLink[dir, shift[dir, coords, -1], 
               getLink[dir, shift[dir, coords, -1]].ConjugateTranspose[u]],
       {dir, nd}];

setRandomGauge::usage = "Perform a random gauge transform.";
setRandomGauge[] :=
 Do[Block[
   {coords = latticeCoordinates[i],
    vt = randomSUMatrix[]},
   If[!boundarySiteQ[coords],
      Do[
          setLink[dd, shift[dd, coords, -1],
	          getLink[dd, shift[dd, coords, -1]].vt];
          setLink[dd, coords, ConjugateTranspose[vt].getLink[dd, coords]],
          {dd, nd}]]],
    {i, latticeVolume[]}];


getBoundaryPhase::usage = "Add up phases going around edge for a given starting point on the boundary.";
getBoundaryPhase::unknown = "Unknown Method";
Options[getBoundaryPhase] = {
    "boundaryMethod" -> "straight", "debug" -> False};
getBoundaryPhase[dir_, coords_, OptionsPattern[]] :=
 Block[{dirt, sign, dt,
        debug = OptionValue["debug"]},
  Scan[(sign = If[#==dir, 1, -1]; If[#!=dir, dirt=#])&,
             transverseBCDirections];
  dt = coords[[dirt]] - (latticeDimensions[[dirt]]+1)/2;
  If[debug, Print["dt=", dt, ", sign=", sign]];
  Which[
      OptionValue["boundaryMethod"] === "discontinuous",
      (* "discontinuous" and "matchBoundary" should be identical
        if the phase on each boundary link is proportional to
        the angle about the middle of the lattice
        subtended by that link. *)
      Block[{fraction = ArcTan[Abs[dt],
                (latticeDimensions[[dir]] - 1)/2]/Pi},
            If[debug, Print[{N[fraction], dt,
                             (latticeDimensions[[dir]] - 1)/2}]];
            DiagonalMatrix[Exp[
                I sign*Sign[dt]*transverseBCPhases*fraction]]],
      OptionValue["boundaryMethod"] === "smooth",
      Block[{fraction = ArcTan[latticeDimensions[[dirt]] -1,
                               latticeDimensions[[dir]] - 1]/Pi},
            fraction *= 2*dt/(latticeDimensions[[dirt]] -1);
            If[debug, Print[{N[fraction], dt,
                             (latticeDimensions[[dir]] - 1)/2}]];
            DiagonalMatrix[Exp[
                I sign*transverseBCPhases*fraction]]],
      OptionValue["boundaryMethod"] === "matchBoundary",
      Block[{y = coords, v = IdentityMatrix[nc]},
       Which[dt > 0,
             While[y[[dirt]] < latticeDimensions[[dirt]],
                   v = v.getLink[dirt, y];
                   y[[dirt]] += 1];
             While[y[[dir]] < latticeDimensions[[dir]],
                   v = v.getLink[dir, y];
                   y[[dir]] += 1];
             While[y[[dirt]] > coords[[dirt]],
                   y[[dirt]] -= 1;
                   v = v.ConjugateTranspose[getLink[dirt, y]]],
             dt < 0,
             While[y[[dirt]] > 1,
                   y[[dirt]] -= 1;
                   v = v.ConjugateTranspose[getLink[dirt, y]]];
             While[y[[dir]] < latticeDimensions[[dir]],
                   v = v.getLink[dir, y];
                   y[[dir]] += 1];
             While[y[[dirt]] < coords[[dirt]],
                   v = v.getLink[dirt, y];
                   y[[dirt]] += 1]];
       v],
      OptionValue["boundaryMethod"] === "straight",
      (* Axial gauge transform without committing any
        crimes at the boundaries. *)
      Block[{y = coords, v = IdentityMatrix[nc]},
       While[y[[dir]] < latticeDimensions[[dir]],
             v = v.getLink[dir, y];
             y[[dir]] += 1];
       v],
      True,
      Message[getBoundaryPhase::unknown];
      Return[$Failed]]]/;
 latticeBC == "TRANS_GAUGEBC" && coords[[dir]] == 1 &&
 MemberQ[transverseBCDirections, dir];


setAxialGauge::usage = "Set axial gauge for the current lattice configuration.
That is, links along any Polyakov loop in direction dir are constant and
diagonal.  In the case \"center\" -> True, set the gauge relative to the
nearest element of the center of the group for each axial-direction link.
Option \"abelian\" -> True rotates the field so that it is diagonal.";
Options[setAxialGauge] = Join[{"center" -> False, "abelian" -> False,
                               "debug" -> False},
                              Options[getBoundaryPhase]];
setAxialGauge::options = "Invalid combination of options.";
setAxialGauge[dir_, opts:OptionsPattern[]] :=
 Block[{popts = Apply[Sequence, FilterRules[{opts},
                      Options[getPhases]]],
        gbopts = Apply[Sequence, FilterRules[{opts},
                       Options[getBoundaryPhase]]],
        debug = OptionValue["debug"]},
  Do[Block[{coords = latticeCoordinates[i]},
   If[coords[[dir]] == 1 && !boundarySiteQ[coords],
    Block[{v = IdentityMatrix[nc], ww, vt, ct, phases, centerValues},
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
        {phases, ww} = Take[getPhases[v, True, popts], 2];
        (* Convert to an SU(N) matrix *)
        ww[[-1]] /= Det[ww];
        v = DiagonalMatrix[Exp[I*phases/latticeDimensions[[dir]]]],
        v = SUPower[v, 1/latticeDimensions[[dir]]]];
     If[debug, Print["Root:", v]];
     Do[
      (* Apply gauge choice to each site *)
      ct = shift[dir, coords, j - 1];
      If[j == 1,
         If[OptionValue["abelian"],
            vt = Transpose[ww],
            Continue[]],
         (* Solving the linear system instead of using the Hermitian
           conjugate decreases the floating point error. *)
         vt = LinearSolve[getLink[dir, shift[dir, ct, -1]], v]];
      If[OptionValue["center"],
         vt *= Exp[I 2 Pi centerValues[[j-1]]/nc]]; 
      Do[setLink[dd, shift[dd, ct, -1],
                 getLink[dd, shift[dd, ct, -1]].vt];
	 setLink[dd, ct,
                 (* Solving the linear system instead of using the
                   Hermitian conjugate decreases the floating point error. *)
                 LinearSolve[vt, getLink[dd,ct]]],
	 {dd, nd}],
      {j, latticeDimensions[[dir]]}];
     (* Verify that it all worked out *)
     If[debug,
	If[False,
           Do[Print[{j, Chop[getLink[dir, shift[dir, coords, j - 1]]]}],
	      {j, latticeDimensions[[dir]]}]];
      Print["Final link diff:",
	    Chop[v - getLink[dir, shift[dir, coords, -1]]]]]]];
   (*
     Special case of fixed boundary fields.
     The site is on a boundary but the link is not.
     Apply axial gauge until the next boundary is reached.
    *)
   If[boundarySiteQ[coords] && !boundaryLinkQ[dir, coords],
    Block[{ui = IdentityMatrix[nc], up, gidagger, gf, gg, vt, y,
           phases, count, ww},
     If[OptionValue["center"] || OptionValue["abelian"],
        Message[setAxialGauge::options];
        Return[]];
     (* Initial values *)
     up = getBoundaryPhase[dir, coords, gbopts];
     y = coords; count = 0;
     While[!boundaryLinkQ[dir, y],
           ui = ui.getLink[dir, y];
           count += 1; y = shift[dir, y]];
     If[debug, Print["Full product:", ui]];
     gidagger = SUPower[up, 0.5].ConjugateTranspose[SUPower[ui, 0.5]];
     gf = ConjugateTranspose[SUPower[ui, 0.5]].SUPower[up, 0.5];
     gg = SUPower[up, 1/count];
     setLink[dir, coords, gidagger.getLink[dir, coords]];
     y = shift[dir, coords];
     (* Apply gauge choice to each interior site *)
     While[
         !boundaryLinkQ[dir, y],
         (* Solving the linear system instead of using the Hermitian
           conjugate decreases the floating point error. *)
         vt = LinearSolve[getLink[dir, shift[dir, y, -1]], gg];
         Do[
             setLink[dd, shift[dd, y, -1],
                     getLink[dd, shift[dd, y, -1]].vt];
	     setLink[dd, y,
                     LinearSolve[vt, getLink[dd, y]]],
	     {dd, nd}];
         y = shift[dir, y]];
     y = shift[dir, y, -1];
     setLink[dir, y, getLink[dir, y].gf];
     (* Verify that it all worked out *)
     If[debug,
        Print["Final link diff:",
	      SUNorm[LinearSolve[gg, getLink[dir, y]]]]]]]],
     {i, latticeVolume[]}]];

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
   (* This instance of ParallelTable responsible for errors
     when calculating observables along trajectories.
     First attempt:  don't specify Method *)
   Do[Scan[update, ParallelTable[
       Block[{coords = latticeCoordinates[i], sum, gauge,
              normSum = 0.0},
             If[Not[ListQ[coords]],
                Messsage[General::list, coords, 0]];
             If[Mod[Total[coords], 2] == cb && !boundarySiteQ[coords],
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
       {i, latticeVolume[]} (*, Method->"CoarsestGrained"*)], {3}], {cb, 0, 1}];
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
           If[coords[[dir1]] == 1 && Mod[Total[coords], 2] == cb &&
              !boundarySiteQ[coords],
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
