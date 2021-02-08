plotActionLimits[] := {0, Min[2 n^2/3, 15]};
plotActionPlane[dir1_, dir2_, anchor_, range_] :=
 Block[{coords = anchor},
  ListPlot3D[
   Flatten[Table[coords[[dir1]] = x1;
     coords[[dir2]] = x2; {x1 + 1/2, x2 + 1/2,
	   beta (1 - Re[plaquette[dir1, dir2, wrapIt[coords]]]/nc)},
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
            Re[plaquette[dir0, If[i == 1, dir1, dir2], wrapIt[coords]]]/nc)},
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

Options[plotPlaquetteCorrelations] = Options[ErrorListPlot];
plotPlaquetteCorrelations[corr_, opts:OptionsPattern[]] := 
 (* Bug in ErrorListPlot:  x-values must be machine numbers *)
 ErrorListPlot[Map[
     {{N[#[[1]]], #[[2]]}, ErrorBar[ #[[3]]]}&,
       (* Order data sets by the value of the last index.
         One could also use GroupBy[] to do this. *)
       Table[Select[corr, Last[#]==i&],
             {i, Union[Map[Last, corr]]}], {2}],
  opts,
  PlotRange -> {{0, All}, All}, Axes -> False, Frame -> True, 
  FrameLabel -> {"separation", "correlation"},
  PlotLegends -> Placed[SwatchLegend[
      (* Data set order is given by the numbering in
        the compiled function "orient" above. *)
      {"same plane", "same orientation, stacked",
       "same orientation, offset", "different orientation"}], {0.7, 0.75}],
  Epilog -> {Dashing[0.01], 
             Line[{{0, 0}, {1/2 + Max[Map[First, corr]], 0}}]}];

Options[plotStringModelFit] = Join[Options[ErrorListPlot],
                                   Options[stringModel],
                                   {"maxx" -> All}];
plotStringModelFit[tallyData_, opts:OptionsPattern[]] :=
 Block[{lowerCutoff = OptionValue["lowerCutoff"], diffs,
        (* Sin[2*theta]^2 *)
        sin2theta2 = Function[{x, y}, (2*x*y/(x^2+y^ 2))^2],
        offAxis = Function[dir, Apply[sin2theta2, Drop[#, {dir}]]&],
        nearest, maxx, miny, maxy, nn, ff, cd = ColorData[1]}, 
  nn = Select[Normal[tallyData], (True || Norm[#[[1, 1]]] > 0) &]; 
  ff = stringModel[tallyData,
       Apply[Sequence, FilterRules[{opts}, Options[stringModel]]],
                   printResult -> True];
  nearest = Length[nn[[1, 1]]] - 1;
  maxx = If[OptionValue["maxx"] === All,
            Max[Map[(Norm[#[[1, 1]]]*#[[1,-1]])&, nn]] +
            Min[latticeDimensions]/2,
            OptionValue["maxx"]];
  maxy = Max[Map[#[[2, 1]] + 1.5 #[[2, 2]] &, nn]]; 
  miny = Min[0, Min[Map[#[[2, 1]] - 1.5 #[[2, 2]] &, nn]]];
  (* Print["limits", {nearest, maxx, miny, maxy, ff}]; *)
  Show[ErrorListPlot[
      Map[
          {{N[Norm[#[[1, 1]]]]*#[[1,-1]], #[[2, 1]]}, ErrorBar[#[[2, 2]]]}&,
          (* Group by ll value *)
          GatherBy[nn, #[[1, -1]]&], {2}],
      Apply[Sequence, FilterRules[{opts}, Options[ErrorListPlot]]],
      PlotRange -> {{0, maxx}, {miny, maxy}}, 
      Epilog -> {Text[
          ColumnForm[
              Join[Map[
                  Row[{#[[1, 1]], " = ", 
                       valueError[#[[1, 2]], #[[2]]]}]&, 
                     Transpose[
                         ff[{"BestFitParameters", 
                             "ParameterErrors"}]]],
                   {Row[{"\[Beta] = ", N[beta], 
                         ", \!\(\*StyleBox[\"N\",\nFontSlant->\"Italic\"]\) = ", nc}],
                    Row[{Row[latticeDimensions, "\[Times]"], "lattice"}, 
                        " "]}]], {maxx, maxy}, {1.1, 1.1}]},
      Prolog -> {{Gray, Line[{{0, 0}, {maxx, 0}}]}, 
                 If[lowerCutoff > 0,
                    {Gray, 
                     Text["cutoff", {lowerCutoff*Min[latticeDimensions], 0},
                          {-1.5, -1}, {0, 1}],
                     Map[
                         Line[{{lowerCutoff*#, miny},
                               {lowerCutoff*#, maxy 2/3}}]&,
                             Union[latticeDimensions]]},
                    Nothing]},
      PlotStyle -> cd,
      Axes -> False, Frame -> True,
      FrameLabel -> {"area enclosed (plaquettes)", 
                     "Polyakov loop correlator"}], 
      Table[
       Block[{face = latticeDimensions, vertices},
        face[[dir0]] = 1;
        vertices = makeVertices[face];
        { (* off-axis case
            Last[vertices] is a highest norm corner. *)
         ParametricPlot[
          Block[{xx = x Last[vertices]/Norm[Last[vertices]]},
                (* Reconstruct the mirrors for this separation.
                  It would be safer to pull this from nn? *)
           diffs = Map[Drop[#, {dir0}]&,
                           Take[SortBy[Map[#-xx&, vertices],
                                       Norm], nearest]];
           {latticeDimensions[[dir0]]*x,
            Apply[ff["Function"],
                  Append[
                      Flatten[diffs],
                      latticeDimensions[[dir0]]]]}],
          {x, Sqrt[nd-1]-1/2, Norm[Last[vertices]]/2+1/2},
          PlotStyle->{cd[dir0], Opacity[0.5], Thickness[0.004],
                      Dashing[0.02]}],
         Table[If[
             dir1 != dir0,
             (* on axis *)
             ParametricPlot[
              Block[{xx = Table[If[i==dir1, x, 0], {i, nd}]},
                    {latticeDimensions[[dir0]]*x,
                Apply[ff["Function"],
                      Append[
                          Flatten[Map[Drop[#, {dir0}]&,
                                          Take[SortBy[Map[#-xx&, vertices],
                                                      Norm], nearest]]],
                          latticeDimensions[[dir0]]]]}],
              {x, 1/2, latticeDimensions[[dir1]]/2+1/2},
              PlotStyle -> {cd[dir0], Opacity[0.5], Thickness[0.004]},
              PlotRange -> All],
             Nothing],
               {dir1, nd}]}],
       {dir0, nd}]]];


Options[plotWilsonModelFit] = Join[Options[Graphics3D],
                                   Options[wilsonModel],
                                   {}];
plotWilsonModelFit[data0_, diffFlag_?BooleanQ, opts___]:=
Block[{bord = 0.025, dims, ff, min, max,
 (* Switch around dimensions to minimize number of plots.
    This does not affect the model fitting. *)
 data = Normal[
     KeyMap[If[#[[3]] > #[[4]], #[[{2, 1, 4, 3}]], #]&, data0]],
 gopts = Apply[Sequence, FilterRules[{opts}, Options[Graphics3D]]]},
 dims = Union[Map[Take[First[#], -2] &, data]]; 
 ff = wilsonModel[data,
      Apply[Sequence, FilterRules[{opts}, Options[wilsonModel]]],
                  printResult -> True]; 
 If[False, Print["Fit residuals: ", ff["StandardizedResiduals"]]]; 
 If[diffFlag,
    max = Max[Map[Block[{y=#[[2, 1]] - Apply[ff["Function"], #[[1]]]},
                        Max[#[[2,2]]+{y,-y}]]&, data]];
    min = -max;
    GraphicsColumn[
        Apply[Function[{l1, l2}, 
         Block[{sdata = 
                Select[data, (#[[1, 3]] == l1 && #[[1, 4]] == l2) &]}, 
          Show[Graphics3D[
              Map[Block[{fff = #[[2, 1]] - Apply[ff["Function"], #[[1]]]},
                        {{Thickness[0.004], Green, 
                          If[Abs[fff] > #[[2, 2]], 
                             Line[{{#[[1, 1]], #[[1, 2]], 0},
                                   {#[[1, 1]], #[[1, 2]], 
                                    If[fff > 0, fff - #[[2, 2]],
                                       fff + #[[2, 2]]]}}], 
                             Nothing]},
                         {Thickness[0.008], Blue, 
                          Line[{{#[[1, 1]], #[[1, 2]], 
                                 fff - #[[2, 2]]},
                                {#[[1, 1]], #[[1, 2]], 
                                 fff + #[[2, 2]]}}]}}] &, sdata],
              gopts,
              Axes -> True, 
              AxesLabel -> {Subscript["w", 1], Subscript["w", 2], "diff"}, 
              BoxRatios -> {1, l2/l1, 0.5}, 
              PlotLabel -> Row[{l1, "\[Times]", l2, " lattice slice"}], 
              PlotRange -> {min, max}], 
  Plot3D[0, {w1, Min[Map[#[[1, 1]] &, sdata]] - bord, 
             Max[Map[#[[1, 1]] &, sdata]] + bord},
         {w2, Min[Map[#[[1, 2]] &, sdata]] - bord, 
          Max[Map[#[[1, 2]] &, sdata]] + bord}, PlotRange -> All, 
         RegionFunction -> (l1 != l2 || #1 >= #2 &), 
         ClippingStyle -> None, PlotStyle -> {Opacity[0.5]}]]]], 
              dims, {1}]],
    min = Min[Map[(#[[2, 1]] - #[[2, 2]] 1.5) &, data]];
    max = Max[Map[(#[[2, 1]] + #[[2, 2]] 1.5) &, data]]; 
    GraphicsColumn[
     Apply[Function[{l1, l2}, 
     Block[{sdata = 
        Select[data, (#[[1, 3]] == l1 && #[[1, 4]] == l2) &]}, 
      Show[Graphics3D[
        Map[{{Thickness[0.004], Gray, 
            Line[{{#[[1, 1]], #[[1, 2]], 
               min}, {#[[1, 1]], #[[1, 
                2]], #[[2, 1]] - #[[2, 2]]}}]}, {Thickness[0.008], 
            Blue, Line[{{#[[1, 1]], #[[1, 
                2]], #[[2, 1]] - #[[2, 2]]}, {#[[1, 1]], #[[1, 
                  2]], #[[2, 1]] + #[[2, 2]]}}]}} &, sdata],
        gopts,
        Axes -> True, 
        AxesLabel -> {Subscript["w", 1], Subscript["w", 2], 
          "\[ScriptCapitalW]"}, 
        BoxRatios -> {1, l2/l1, 0.5}, FaceGrids -> {{0, 0, -1}}, 
        PlotLabel -> Row[{l1, "\[Times]", l2, " lattice slice"}], 
        PlotRange -> {min, max}], 
        Plot3D[
            ff["Function"][w1, w2, l1, l2],
            {w1, Min[Map[#[[1, 1]] &, sdata]] - bord, 
             Max[Map[#[[1, 1]] &, sdata]] + bord},
            {w2, Min[Map[#[[1, 2]] &, sdata]] - bord, 
             Max[Map[#[[1, 2]] &, sdata]] + bord},
            PlotRange -> All, 
            RegionFunction -> (l1 != l2 || #1 >= #2 &), 
            ClippingStyle -> None,
            PlotStyle -> If[True, None, Directive[Opacity[0.4]]]]]]], 
      dims, {1}]]]];

plotDist::usage = "Distribution of complex numbers in a unit circle in the complex plane, showing Im[z]<0 part, and the center of the gauge group.";
Options[plotDist] := Options[Histogram3D];
plotDist[data_, op_:"1", all_:False, opts:OptionsPattern[]] := (
  Print["{mean, Re stdev, Im stdev, count}: ",
        {Mean[Re[data]], StandardDeviation[Re[data]],
         Chop[Sqrt[Mean[Map[Im[#]^2&, data]]]], Length[data]}];
  Show[{Histogram3D[Map[{Re[#], -Im[#]} &, data], {0.05},
     opts,
     PlotRange -> {{-1, 1}, {-1, If[all, 1, 0]}, {0, All}}, 
                    BoxRatios -> {1, If[all, 1, 0.5], 1/3}, 
     AxesLabel -> {"Re Tr(\[ScriptCapitalO])",
                   "Im Tr(\[ScriptCapitalO])", None}, 
     FaceGrids -> If[all, None, {{0, 1, 0}}]], 
     Graphics3D[{(*{Red, Thickness[0.01], 
       Table[Line[{{0, 0, 0.1}, 
                   {Cos[2 Pi k/nc], -Sin[2 Pi k/nc], 1.5}}], {k, 0, nc/2}]},*)
         Table[
             If[all || Sin[2 Pi k/nc]<=0,
                Text[z^k, {Cos[2 Pi k/nc], Sin[2 Pi k/nc], 0}, {1-2 k/nc, -2}],
                Nothing],
             {k, 0, nc-1}],
         {Blue, Thickness[0.005],
          Table[Line[Map[{Re[#], Im[#], 1.5}&,
            Table[stringOperator[
                centerSUMatrix[k1].SUPower[centerSUMatrix[k2], x], op],
                  {x, 0, 1, 1/(10 nc)}]]], {k2, 1, nc-1}, {k1, 0, k2}]},
         {Hue[0.5, .3], 
          Cylinder[{{0, 0, 0.05}, {0, 0, -0.1}}, 1]}}]}])/;op != "phases";

Options[plotPhases] := Join[{"centerLabels" -> True, "bins" -> {Pi/20}},
                            If[nc==2, Options[Histogram],
                               Options[Histogram3D]]];
plotPhases[data_, all_:False, opts:OptionsPattern[]] :=
    Block[{xyEdge = 0.2, dy = Sqrt[Length[data]]/100, z, \[Lambda], basis,
           reflect = If[#1, -Abs[#2], #2]&,
           conj = Function[k, If[k==0, 0, nc-k], {Listable}],
        centerPhases = Function[k, Sort[
            cleanPhases[Table[2 Pi k/nc, {nc}]], Greater], {Listable}]},
    Print["{{mean, stdev}, ...}: ",
          Transpose[{Mean[data], StandardDeviation[data]}]];
    Which[
        nc == 2,
        Histogram[Map[First, data],
                  Apply[Sequence, FilterRules[{opts}, Options[Histogram]]],
                  PlotRange-> {{-Pi, Pi}, {0, All}},
                  Axes -> False, Frame -> True,
                  FrameLabel -> {Subscript[\[Lambda], 1], "count"}
                  ],
        nc == -3,
        Block[{
            grid = Table[x, {x, -2 Pi, 2 Pi, Pi/nc}],
            ticks1 = Table[
                If[Mod[x, Pi] == 0, {x, x, 0.03}, {x, "", 0.015}],
                {x, -2 Pi, 2 Pi, Pi/nc}],
            ticks2 = Table[{x, x, 0.03}, {x, -2 Pi, 2 Pi, Pi/nc}]},
          Show[Histogram3D[
              Map[{#[[1]], #[[2]]} &, data], {Pi/30},
              Apply[Sequence, FilterRules[{opts}, Options[Histogram3D]]],
              PlotRange -> {{0-xyEdge, 4 Pi/3+xyEdge},
                            {-2 Pi/3-xyEdge, 2 Pi/3+xyEdge}, {0, All}}, 
              BoxRatios -> {1, 1, 1/3},
              Ticks -> {ticks1, ticks2, Automatic},
              AxesLabel -> {Subscript[\[Lambda], 1],
                            Subscript[\[Lambda], 2], None}, 
              FaceGrids -> {{{0, 0, -1}, {grid, grid}}}], 
               Graphics3D[{
                   {Text[z^2, {4 Pi/3, -2 Pi/3, 0}, {-1, 1}],
                    Text[z^1, {2 Pi/3, 2 Pi/3, 0}, {-1, -1}],
                    Text[z^0, {0, 0, 0}, {1, -1}]},
                   {Blue, Thickness[0.005],
                    Line[{{0, 0, dy}, {4 Pi/3, -2 Pi/3, dy},
                          {2 Pi/3, 2 Pi/3, dy}, {0, 0, dy}}]}}]]],
        True,
        basis = Table[Subscript[\[Lambda], i], {i, nc}];
        GraphicsColumn[Flatten[Table[
            Block[{v0, v1, v2, b1, b2, dd, oo, perm,
                   max1, min2, max2, grid1, grid2, ticks1, ticks2},
                  oo = Order[{k0, k1, k2}, Sort[conj[{k0, k1, k2}]]];
               If[oo < 0,
                  Print["Skipping ", {k0, k1, k2}]; Nothing,
                  perm = Which[
                      oo==0 && k2 == conj[k2], {k2, k0, k1},
                      oo==0 && k1 == conj[k1], {k1, k0, k2},
                      True, {k0, k1, k2}];
                  {v0, v1, v2} = centerPhases[perm];
                  v1 -= v0; v2 -= v0;
                  b1 = Simplify[(v1 + v2)/(2 Pi)];
                  b2 = Simplify[(v1 - v2)/(2 Pi)];
                  dd = Map[{(#-v0).b1, reflect[oo==0 && !all, (#-v0).b2]}&,
                           data];
                  If[True,
                     Print[Chop[{{k0, k1, k2}, v0, {v1, v2}, {b1, b2}}]]];
                  max1 = Max[v1.b1, v2.b1];
                  min2 = Min[v1.b2, v2.b2];
                  max2 = If[oo==0 && !all, 0, Max[v1.b2, v2.b2]];
                  grid1 = Table[x, {x, 0, max1, Pi/nc}];
                  grid2 = Table[x, {x, min2, max2, Pi/nc}];
                  ticks1 = Table[
                      If[Mod[x, Pi] == 0, {x, x, 0.03}, {x, "", 0.015}],
                      {x, 0, max1, Pi/nc}];
                  ticks2 = Table[
                      If[Mod[x, Min[Pi, Abs[min2]]] == 0, {x, x, 0.03},
                         {x, "", 0.015}], {x, min2, max2, Pi/nc}];
                  Show[Histogram3D[
                      dd, OptionValue["bins"],
                      Apply[Sequence,
                            FilterRules[{opts}, Options[Histogram3D]]],
                      PlotRange -> {{0, max1+xyEdge},
                                    {min2-xyEdge,
                                     If[oo == 0 && !all, 0, max2+xyEdge]},
                                    {0, All}},
                      BoxRatios -> {1, If[oo==0 && !all, 1/2, 1]*2/Sqrt[3],
                                    1/2},
                      Ticks -> {ticks1, ticks2, {10000, 20000}},
                      FaceGrids -> Join[
                          If[oo==0 && !all,
                             {{{0, 1, 0}, {grid1, Automatic}}}, {}],
                          {{{0, 0, -1}, {grid1, grid2}}}],
                      AxesLabel -> If[True,
                                      {Simplify[
                                          (b1-Table[1, {nc}] Last[b1]).basis],
                                       Simplify[
                                           (b2-Table[1, {nc}] Last[b2]).basis],
                                       None},
                                      None],
                      AspectRatio -> 1/GoldenRatio], 
                    Graphics3D[{
                       (* May need to turn off explicitly if plot range
                         is modified *)
                       If[OptionValue["centerLabels"],
                          {If[oo == 0 && !all, Nothing,
                              Text["z"^perm[[2]], {v1.b1, v1.b2, 0}, {1.5, -2}]],
                           Text["z"^perm[[3]], {v2.b1, v2.b2, 0}, {0.5, -2}],
                           Text["z"^perm[[1]], {0, 0, 0}, {1.5, -1}]}, Nothing],
                       {Blue, Thickness[0.005],
                        If[oo==0 && !all,
                           Line[{{(v1+v2).b1/2, (v1+v2).b2/2, dy},
                                 {v2.b1, v2.b2, dy}, {0, 0, dy}}],
                           Line[{{0, 0, dy}, {v1.b1, v1.b2, dy},
                                 {v2.b1, v2.b2, dy}, {0, 0, dy}}]]}}]]]],
            {k0, 0, nc-3}, {k1, k0+1, nc-2}, {k2, k1+1, nc-1}]]]
    ]];

plotPhaseFunction::usage = "Function version of plotPhases[], doesn't do anything useful for nc>3.";
Options[plotPhaseFunction] := If[nc==2, Options[Plot],
                          Options[Plot3D]];
plotPhaseFunction[f_, all_:False, opts:OptionsPattern[]] :=
    Block[{xyEdge = 0.2, dy = 0.02, z, \[Lambda], basis,
           reflect = If[#1, -Abs[#2], #2]&,
           conj = Function[k, If[k==0, 0, nc-k], {Listable}],
        centerPhases = Function[k, Sort[
            cleanPhases[Table[2 Pi k/nc, {nc}]], Greater], {Listable}],
           grid = Table[x, {x, -2 Pi, 2 Pi, Pi/nc}],
           ticks1 = Table[
               If[Mod[x, Pi] == 0, {x, x, 0.03}, {x, "", 0.015}],
               {x, -2 Pi, 2 Pi, Pi/nc}],
           ticks2 = Table[{x, x, 0.03}, {x, -2 Pi, 2 Pi, Pi/nc}]},
    Which[
        nc == 2,
        Plot[f[{x, -x}], {x, 0, Pi},
                  opts,
                  PlotRange-> {0, All},
                  Axes -> False, Frame -> True,
                  FrameLabel -> {Subscript[\[Lambda], 1], "probability"}
                  ],
        nc == -3,
        Show[Plot3D[
            f[{x, y, -x-y}], {x, 0, 4 Pi/3}, {y, 0, 4 Pi/3},
            opts,
            PlotRange -> {{0-xyEdge, 4 Pi/3+xyEdge},
                          {-2 Pi/3-xyEdge, 2 Pi/3+xyEdge}, {0, All}}, 
            BoxRatios -> {1, 1, 1/3},
            Ticks -> {ticks1, ticks1, Automatic},
            AxesLabel -> {Subscript[\[Lambda], 1],
                          Subscript[\[Lambda], 2], "probability"}, 
            FaceGrids -> {{{0, 0, -1}, {grid, grid}}}], 
              Graphics3D[{
                  {Text[z^2, {4 Pi/3, -2 Pi/3, 0}, {-1, 1}],
                   Text[z^1, {2 Pi/3, 2 Pi/3, 0}, {-1, -1}],
                   Text[z^0, {0, 0, 0}, {1, -1}]},
                  {Blue, Thickness[0.005],
                   Line[{{0, 0, dy}, {4 Pi/3, -2 Pi/3, dy},
                         {2 Pi/3, 2 Pi/3, dy}, {0, 0, dy}}]}}]],
        True,
        basis = Table[Subscript[\[Lambda], i], {i, nc}];
        GraphicsColumn[Flatten[Table[
            Block[{v0, v1, v2, b1, b2, oo, perm,
                   max1, min2, max2, grid1, grid2, ticks1, ticks2},
                  oo = Order[{k0, k1, k2}, Sort[conj[{k0, k1, k2}]]];
                  If[oo < 0,
                     Print["Skipping ", {k0, k1, k2}]; Nothing,
                     perm = Which[
                         oo==0 && k2 == conj[k2], {k2, k0, k1},
                         oo==0 && k1 == conj[k1], {k1, k0, k2},
                         True, {k0, k1, k2}];
                  {v0, v1, v2} = centerPhases[perm];
                  v1 -= v0; v2 -= v0;
                  b1 = Simplify[(v1 + v2)/(2 Pi)];
                  b2 = Simplify[(v1 - v2)/(2 Pi)];
                  If[True,
                     Print[Chop[{{k0, k1, k2}, oo==0, v0,
                                 {v1, v2}, {b1, b2}}]]];
                  max1 = Max[v1.b1, v2.b1];
                  min2 = Min[v1.b2, v2.b2];
                  max2 = If[oo==0 && !all, 0, Max[v1.b2, v2.b2]];
                  grid1 = Table[x, {x, 0, max1, Pi/nc}];
                  grid2 = Table[x, {x, min2, max2, Pi/nc}];
                  ticks1 = Table[
                      If[Mod[x, Pi] == 0, {x, x, 0.03}, {x, "", 0.015}],
                      {x, 0, max1, Pi/nc}];
                  ticks2 = Table[
                      If[Mod[x, Min[Pi, Abs[min2]]] == 0, {x, x, 0.03},
                         {x, "", 0.015}], {x, min2, max2, Pi/nc}];
                  Show[Plot3D[
                      f[v0 + x b1/b1.b1 + y b2/b2.b2],
                      {x, 0, max1 + xyEdge},
                      {y, min2 - xyEdge,
                       If[oo == 0 && !all, 0, max2 + xyEdge]},
                      opts,
                      PlotRange -> {0, All},
                      BoxRatios -> {1, If[oo==0 && !all, 1/2, 1]*2/Sqrt[3],
                                    1/2},
                      Ticks -> {ticks1, ticks2, Automatic},
                      RegionFunction -> Function[{x, y},
                             Block[{vv = v0 + x b1/b1.b1 + y b2/b2.b2},
                               vv[[1]]>vv[[2]] &&
                               vv[[2]]>vv[[3]] && vv[[1]]-vv[[3]]<2 Pi]],
                      AxesLabel -> If[True,
                                      {Simplify[
                                          (b1-Table[1, {nc}] Last[b1]).basis],
                                       Simplify[
                                           (b2-Table[1, {nc}] Last[b2]).basis],
                                       None},
                                      None],
                      AspectRatio -> 1/GoldenRatio], 
                   Graphics3D[{
                       {If[oo == 0 && !all, Nothing,
                           Text["z"^perm[[2]], {v1.b1, v1.b2, 0},
                                {1.5, -2}]],
                        Text["z"^perm[[3]], {v2.b1, v2.b2, 0},
                             {0.5, -2}],
                        Text["z"^perm[[1]], {0, 0, 0}, {1.5, -1}]},
                       {Blue, Thickness[0.005],
                        If[oo==0 && !all,
                           Line[{{(v1+v2).b1/2, (v1+v2).b2, dy},
                                 {v2.b1, v2.b2, dy}, {0, 0, dy}}],
                           Line[{{0, 0, dy}, {v1.b1, v1.b2, dy},
                                 {v2.b1, v2.b2, dy}, {0, 0, dy}}]]}}],
                       (* No idea why, but FaceGrids inside Plot3D[]
                         are ignored. *)
                       FaceGrids -> Join[
                           If[oo==0 && !all,
                              {{{0, 1, 0}, {grid1, Automatic}}}, {}],
                           {{{0, 0, -1}, {grid1, grid2}}}]]]],
            {k0, 0, nc-3}, {k1, k0+1, nc-2}, {k2, k1+1, nc-1}]]]
    ]];


(*
           Plot gauge field configuration
 *)

getMaxField2[lgf_, stdev_:2] :=
    Map[Min[Max[#], Mean[#] + stdev*StandardDeviation[#]]&,
           {Flatten[Map[Norm[#, "Frobenius"]&, lgf, {2}]],
            Abs[Flatten[lgf]]}];
plotA[a_, center_, maxSize_List] := 
    Append[plotA[a, 0, maxSize], 
           Text[center, {1/2, 0}, {0, 0}]]/;center!=0; 
plotA[a_, 0, {maxSizeF_, maxSize1_}] := {{GrayLevel[0.8],
 Block[{z = Min[1, Norm[a, "Frobenius"]/maxSizeF]/4},
  Polygon[{{0, 0}, {1/2, z}, {1, 0}, {1/2, -z}}]]}, 
  Table[Block[{size = Min[1, Sqrt[Abs[a[[i, j]]]/maxSize1]]/(4 nc),
               x = 1/4 + (i - 1/2)/(2 nc), 
               y = 1/4 - (j - 1/2)/(2 nc)},
    If[nc Abs[a[[i, j]]] >= maxSize, big2++, small2++];
    {Hue[(Arg[a[[i, j]]] +  Pi)/(2 Pi)],(*Print[N[{x,y,size}]];*)
     Rectangle[{x - size, y - size}, {x + size, y + size}]}],
       {i, nc}, {j, nc}]};

plotLinksPlane[lgf_, centers_, dir1_, dir2_, anchor_, range_, size_] := 
    Graphics[{
        Block[{coords = anchor, k,
               value = If[#1 === None, 0, Part[##]]&},
   Flatten[Table[
      coords[[dir1]] = x1;
      coords[[dir2]] = x2;
      k = linearSiteIndex[wrapIt[coords]];
      {If[x1 > 0 && x2 > 0, Point[{x1, x2}], Nothing], 
       Translate[
           plotA[lgf[[dir1, k]], value[centers, dir1, k], range],
           {x1, x2}], 
       Translate[
           Rotate[plotA[lgf[[dir2, k]], value[centers, dir2, k], range], 
                  Pi/2, {0, 0}], {x1, x2}]},
      {x1, 0, latticeDimensions[[dir1]]},
      {x2, 0, latticeDimensions[[dir2]]}], 2]],
            {PointSize[Large], Point[anchor[[{dir1, dir2}]]]}},
          PlotRange -> {{1/2, latticeDimensions[[dir1]] + 1/2},
                        {1/2, latticeDimensions[[dir2]] + 1/2}},
          ImageSize -> size]/;dir1!=dir2;

plotLinksSides[lgf_, centers_, dir0_, dir1_, dir2_, anchor_, range_, size_] := 
 Block[{coords = anchor, k, value = If[#1 === None, 0, Part[##]]&}, 
  Graphics[Flatten[Table[
      coords[[dir1]] = x1;
      coords[[dir2]] = x2;
      k = linearSiteIndex[coords];
      {Point[{x1, x2}], 
       Translate[
           Rotate[plotA[lgf[[dir0, k]], value[centers, dir0, k], range], 
                  Pi/4, {0, 0}], {x1, x2}]},
      {x1, latticeDimensions[[dir1]]}, {x2, latticeDimensions[[dir2]]}], 
     2], PlotRange -> {{1/2, latticeDimensions[[dir1]] + 1/2}, {1/2, 
   latticeDimensions[[dir2]] + 1/2}}, ImageSize -> size]]/;
  dir1!=dir2&&dir0!=dir1&&dir0!=dir2;

Options[plotLinks] = {anchor :> Table[1, {nd}], "stDev" -> 2,
   ImageSize -> All, "centerFlag" -> True};
plotLinks[dir0_, dir1_, dir2_, OptionsPattern[]] := 
 Block[{coords = OptionValue[anchor], max,
        lgf = Map[Chop[-I*SULog[#, "center"-> OptionValue["centerFlag"]]]&,
                      gaugeField, {2}],
        centers = Map[SUNorm[#, "center" -> OptionValue["centerFlag"]][[2]]&,
                            gaugeField, {2}]},
  max = getMaxField2[lgf, OptionValue["stDev"]];
  ListAnimate[Flatten[ParallelTable[
       coords[[dir0]] = x0;
       {plotLinksPlane[lgf, centers, dir1, dir2, coords,
                       max, OptionValue[ImageSize]], 
        plotLinksSides[lgf, centers, dir0, dir1, dir2, coords, 
                       max, OptionValue[ImageSize]]},
       {x0, latticeDimensions[[dir0]]}], 1],
               AnimationRate -> 0.2, AnimationRunning -> False, 
               Deployed -> True]];
