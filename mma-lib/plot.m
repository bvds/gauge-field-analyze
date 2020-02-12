(* See https://mathematica.stackexchange.com/questions/132454 *)
printNumberError[a_, b_] :=
    Block[{exponent = Floor[Max[RealExponent[{a, b}]]],
   sas =  If[Floor[RealExponent[#1]] >= Floor[RealExponent[#2]] - 1,
             SetAccuracy[#1, Accuracy[SetPrecision[#2, 2]]], 0]&},
  If[exponent < 6 && exponent > -4,
     NumberForm[Row[{sas[a, b], "±", SetPrecision[b, 2]}, " "],
                ExponentFunction -> (Null &)],
     Block[{aa = a/10^exponent, bb = b/10^exponent},
           Row[{"(", sas[aa, bb], "±", SetPrecision[bb, 2], ")",
                "×", Superscript[10, exponent]}, " "]]]];
Format[valueError[a_, b_]] := printNumberError[a, b];

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
                                   Options[stringModel]];
plotStringModelFit[tallyData_, opts:OptionsPattern[]] :=
 (Format[a2sigma] = Row[{Style[a, Italic]^2, \[Sigma]}];
  Block[{lowerCutoff = OptionValue["lowerCutoff"],
         nearest, maxx, miny, maxy, nn, ff}, 
    nn = Select[Normal[tallyData], #[[1, 1]] > 0 &]; 
    ff = stringModel[tallyData,
      Apply[Sequence, FilterRules[{opts}, Options[stringModel]]],
                     printResult -> True];
    nearest = Length[nn[[1,1]]] - 1;
    maxx = Max[Map[#[[1, 1]]&, nn]] + Min[latticeDimensions]/2;
    maxy = Max[Map[#[[2, 1]] + 1.5 #[[2, 2]] &, nn]]; 
    miny = Min[0, Min[Map[#[[2, 1]] - 1.5 #[[2, 2]] &, nn]]]; 
    Show[ErrorListPlot[
     Map[{{N[#[[1, 1]]], #[[2, 1]]}, ErrorBar[#[[2, 2]]]} &, nn],
     Apply[Sequence, FilterRules[{opts}, Options[ErrorListPlot]]],
     PlotRange -> {{0, maxx}, {miny, maxy}}, 
     Epilog -> {Text[
         ColumnForm[
             Join[Map[
                 Row[{#[[1, 1]], " = ", 
                      printNumberError[#[[1, 2]], #[[2]]]}] &, 
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
                    Text["cutoff", {lowerCutoff, 0}, {-1.5, -1}, {0, 1}],
                    Line[{{lowerCutoff, miny}, {lowerCutoff, maxy/2}}]},
                   Nothing]},
     Axes -> False, Frame -> True, 
     FrameLabel -> {"area enclosed (plaquettes)", 
                    "Polyakov loop correlator"}], 
     Table[
      Block[{face = latticeDimensions, vertices},
       face[[dir0]] = 1;
       vertices = makeVertices[face];
       {ParametricPlot[
         Block[{xx = x Last[vertices]/Norm[Last[vertices]]},
                 {latticeDimensions[[dir0]]*x,
           Apply[ff,
                 Append[
                     Take[Sort[Map[Norm[#-xx]&, vertices], Less], nearest]*
                     latticeDimensions[[dir0]],
                     latticeDimensions[[dir0]]]]}],
         {x, Sqrt[nd-1]-1/2, Norm[Last[vertices]]/2+1/2},
         PlotRange -> All, PlotStyle->{Orange, Dashing[0.02]}],
        Table[If[
            dir1 != dir0,
            ParametricPlot[
             Block[{xx = Table[If[i==dir1, x, 0], {i, nd}]},
               {latticeDimensions[[dir0]]*x,
              Apply[ff,
                    Append[
                        Take[Sort[Map[Norm[#-xx]&, vertices], Less], nearest]*
                        latticeDimensions[[dir0]],
                        latticeDimensions[[dir0]]]]}],
             {x, 1/2, latticeDimensions[[dir1]]/2+1/2},
             PlotRange -> All],
            Nothing],
              {dir1, nd}]}],
      {dir0, nd}]]]);

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

Options[plotPhases] := If[nc==2, Options[Histogram],
                          Options[Histogram3D]];
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
                  opts,
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
              opts,
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
                      dd, {Pi/20},
                      opts,
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
                       {If[oo == 0 && !all, Nothing,
                           Text["z"^perm[[2]], {v1.b1, v1.b2, 0}, {1.5, -2}]],
                           Text["z"^perm[[3]], {v2.b1, v2.b2, 0}, {0.5, -2}],
                        Text["z"^perm[[1]], {0, 0, 0}, {1.5, -1}]},
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

plotLink[u_, maxSize_] := 
 Block[{a = Chop[-I SULog[u, "center" -> True]], 
   norm = SUNorm[u, "center" -> True]}, 
  Append[plotA[a, maxSize], 
         If[norm[[2]] != 0, Text[norm[[2]], {1/2, 0}, {0, 0}], Nothing]]]; 
plotA[a_, maxSize_] := {{GrayLevel[0.8],
  Block[{z = Min[1, Sqrt[2]*Norm[Flatten[a]]/maxSize]/4}, 
    Polygon[{{0, 0}, {1/2, z}, {1, 0}, {1/2, -z}}]]}, 
  Table[Block[{size = 
      Min[1, Sqrt[Sqrt[2] nc Abs[a[[i, j]]]/maxSize]]/(4 nc), 
     x = 1/4 + (i - 1/2)/(2 nc), 
     y = 1/4 - (j - 1/2)/(2 nc)}, {Hue[(Arg[a[[i, j]]] + 
         Pi)/(2 Pi)],(*Print[N[{x,y,size}]];*)
     Rectangle[{x - size, y - size}, {x + size, y + size}]}],
        {i, nc}, {j, nc}]};

plotLinksPlane[dir1_, dir2_, anchor_, range_, size_, groupFlag_] := 
  Block[{coords = anchor, 
         plotIt = If[groupFlag, plotLink, plotA]}, 
   Graphics[Flatten[Table[coords[[dir1]] = x1;
      coords[[dir2]] = 
       x2; {If[x1 > 0 && x2 > 0, Point[{x1, x2}], Nothing], 
       Translate[
        plotIt[getLink[dir1, wrapIt[coords]], range], {x1, x2}], 
       Translate[
        Rotate[plotIt[getLink[dir2, wrapIt[coords]], range], 
         Pi/2, {0, 0}], {x1, x2}]}, {x1, 0, 
       latticeDimensions[[dir1]]}, {x2, 0, 
       latticeDimensions[[dir2]]}], 2], 
    PlotRange -> {{1/2, latticeDimensions[[dir1]] + 1/2}, {1/2, 
       latticeDimensions[[dir2]] + 1/2}}, ImageSize -> size]]/;dir1!=dir2;
plotLinksSides[dir0_, dir1_, dir2_, anchor_, range_, size_, 
   groupFlag_] := 
  Block[{coords = anchor, 
         plotIt = If[groupFlag, plotLink, plotA]}, 
   Graphics[Flatten[Table[coords[[dir1]] = x1;
      coords[[dir2]] = x2; {Point[{x1, x2}], 
       Translate[
        Rotate[plotIt[getLink[dir0, coords], range], 
         Pi/4, {0, 0}], {x1, x2}]}, {x1, 
       latticeDimensions[[dir1]]}, {x2, latticeDimensions[[dir2]]}], 
     2], PlotRange -> {{1/2, latticeDimensions[[dir1]] + 1/2}, {1/2, 
   latticeDimensions[[dir2]] + 1/2}}, ImageSize -> size]]/;
  dir1!=dir2&&dir0!=dir1&&dir0!=dir2;
Options[plotLinks] = {anchor :> Table[1, {nd}], PlotRange :> 1, 
   ImageSize -> All, groupFlag -> True};
plotLinks[dir0_, dir1_, dir2_, OptionsPattern[]] := 
    Block[{coords = OptionValue[anchor]}, 
   ListAnimate[Flatten[Table[coords[[dir0]] = x0;
      {plotLinksPlane[dir1, dir2, coords, OptionValue[PlotRange], 
        OptionValue[ImageSize], OptionValue[groupFlag]], 
       plotLinksSides[dir0, dir1, dir2, coords, 
        OptionValue[PlotRange], OptionValue[ImageSize], 
        OptionValue[groupFlag]]}, {x0, latticeDimensions[[dir0]]}], 
     1], AnimationRate -> 0.2, AnimationRunning -> False, 
               Deployed -> True]];

