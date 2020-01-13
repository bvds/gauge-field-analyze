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

plotDist::usage = "Distribution of complex numbers in a unit circle in the complex plane, showing Im[z]<0 part, and the center of the gauge group.";
plotDist[data_] := 
  Show[{Histogram3D[Map[{Re[#], -Im[#]} &, data], 
     PlotRange -> {{-1, 1}, {-1, 0.01}, {0, All}}, 
     BoxRatios -> {1, 0.5, 0.25}, 
     AxesLabel -> {"Re(\[CapitalPhi])", "Im(\[CapitalPhi])", None}, 
     FaceGrids -> {{0, 1, 0}}, BaseStyle -> journalStyle], 
    Graphics3D[{{Red, Thickness[0.01], 
       Table[Line[{{0, 0, 0.1}, 
          2 {Cos[2 Pi k/nc], Sin[2 Pi k/nc], 0.1}}], {k, 0, 
         nc - 1}]}, {Hue[0.5, .3], 
       Cylinder[{{0, 0, 0.05}, {0, 0, -0.1}}, 1]}}]}];


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
       latticeDimensions[[dir2]] + 1/2}}, ImageSize -> size]];
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
       latticeDimensions[[dir2]] + 1/2}}, ImageSize -> size]];
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

