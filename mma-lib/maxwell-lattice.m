(* General utilities *)

id[x_] := IdentityMatrix[3][[x]];
(* uses globally defined association "index" *)
addTo[k_?IntegerQ, co_?ListQ, dx_: 1] := 
  If[KeyExistsQ[index, co], aa[[k, index[co]]] += dx];


(* Poisson equation:  2 dimensions *)

makeCoords2[r_] := (maxx = 0; basis = 0; radius = r; 
  coords = Association[];
  index = Association[]; 
  Do[If[Norm[{i, j}] < r, basis += 1; If[maxx < i, maxx = i]; 
                          coords[basis] = {i, j}; index[{i, j}] = basis],
     {i, 0, r}, {j, 0, i}])

makeMatrix2[] := (aa = 
   Block[{x, y}, 
    SparseArray[
     Flatten[Table[{x, y} = coords[i]; 
       Append[Table[
         If[KeyExistsQ[
           index, {x, y} + delta], {i, index[{x, y} + delta]} -> -1, 
          Nothing], {delta, {{-1, 0}, {1, 0}, {0, -1}, {0, 1}}}], {i, 
          i} -> 4], {i, basis}]], basis]]; 
  aa[[index[{0, 0}], index[{1, 0}]]] -= 2; 
  Block[{x, y, epsilon}, 
   Do[{x, y} = coords[k]; 
    If[! KeyExistsQ[index, {x + 1, y}], 
     epsilon = N[Sqrt[radius^2 - y^2] - x]; 
     If[False, Print[{{x, y}, epsilon}]]; aa[[k, k]] += 1/epsilon - 1;
      If[! KeyExistsQ[index, {x, y + 1}], 
      epsilon = N[Sqrt[radius^2 - x^2] - y]; 
      If[False, Print[{{x, y}, epsilon}]]; 
      aa[[k, k]] += 1/epsilon - 1]], {k, basis}]]; 
  Do[If[KeyExistsQ[index, {k, 1}], 
    aa[[index[{k, 0}], index[{k, 1}]]] -= 1]; 
   If[KeyExistsQ[index, {k, k}] && k > 0, 
    aa[[index[{k, k}], index[{k, k - 1}]]] -= 1]; 
   If[KeyExistsQ[index, {k + 1, k}], 
    aa[[index[{k, k}], index[{k + 1, k}]]] -= 1], {k, 0, maxx}]; 
  b = Table[0, {basis}]; b[[index[{0, 0}]]] = 1.0;)

energy2[] := 
 Block[{ee = 0, 
   getvv = (If[KeyExistsQ[index, #], vv[[index[#]]], 0] &)}, 
  Do[ee += If[j == 0, 0.5, 1]*(getvv[{i, j}] - getvv[{i + 1, j}])^2; 
     If[j < i, ee += (getvv[{i, j}] - getvv[{i, j + 1}])^2],
     {i, 0, maxx}, {j, 0, i}]; ee]


(* Poisson equation: 3 dimensions *)

(* Note:  This cannot be used for the 2+1 dimension Wilson loops. *)
makeCoords3[r_] := (maxx = 0; basis = 0; radius = r; 
   coords = Association[];
   index = Association[]; 
   Do[If[Norm[{i, j, k}] < r,
         basis += 1; If[maxx < i, maxx = i]; 
         coords[basis] = {i, j, k}; index[{i, j, k}] = basis],
      {i, 0, r}, {j, 0, i}, {k, 0, j}]);

makeMatrix3[] := (aa = 
   Block[{x, y, z}, 
    SparseArray[
     Flatten[Table[{x, y, z} = coords[i]; 
       Append[Table[
         If[KeyExistsQ[
           index, {x, y, z} + delta], {i, 
            index[{x, y, z} + delta]} -> -1, 
          Nothing], {delta, {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1,
             0}, {0, 0, -1}, {0, 0, 1}}}], {i, i} -> 6], {i, basis}]],
      basis]]; 
  Block[{x, y, z, epsilon}, 
   Do[{x, y, z} = coords[k];
    (* radius border *)
    If[! KeyExistsQ[index, {x + 1, y, z}], 
     epsilon = N[Sqrt[radius^2 - y^2 - z^2] - x]; 
     aa[[k, k]] += 1/epsilon - 1; 
     If[! KeyExistsQ[index, {x, y + 1, z}], 
      epsilon = N[Sqrt[radius^2 - x^2 - z^2] - y]; 
      aa[[k, k]] += 1/epsilon - 1]; 
     If[! KeyExistsQ[index, {x, y, z + 1}], 
      epsilon = N[Sqrt[radius^2 - x^2 - y^2] - z]; 
      aa[[k, k]] += 1/epsilon - 1]];
    (* top face *)
    If[y == z, addTo[k, {x, y, z - 1}]; addTo[k, {x, y + 1, z}]];
    (* bottom face *)
    If[z == 0 , addTo[k, {x, y, 1}]];
    (* side face *)  
    If[x == y, addTo[k, {x + 1, y, z}]; addTo[k, {x, y - 1, z}]];
    (* Edge *)
    If[y == z == 0, addTo[k, {x, 1, 0}, 2]]; 
    If[z == x == y, addTo[k, {x, y, z - 1}]; addTo[k, {x + 1, y, z}]];
    If[x == 0 && y == 0 && z == 0, addTo[k, {1, 0, 0}, 3]],
      {k, basis}]];
  b = Table[0, {basis}];
  b[[index[{0, 0, 0}]]] = 1.0;)


(*  Magnetostatics:  3 dimensions, no symmetries *)

makeCoordsW[r_, gauge_: "Coulomb"] := (iDir = 3; basis = 0; 
  radius = r; coords = Association[];
  index = Association[];
  Do[Block[{x = {i, j, k}, dx}, dx = x + id[dir]/2 - id[iDir]/2;
    (* Exclude links with midpoints on the boundary *)
    If[Norm[dx] < radius && Not[IntegerQ[gauge] && dir == gauge], 
     basis += 1; coords[basis] = {dir, x}; 
     index[{dir, x}] = basis]], {dir, 3}, {i, -Floor[radius], 
    Ceiling[radius]}, {j, -Floor[radius], 
    Ceiling[radius]}, {k, -Floor[radius], Ceiling[radius]}]; 
  gBasis = basis;
  If[gauge === "Coulomb",
     Do[Block[{x = {i, j, k}},
     (* Apply Coulomb gauge constraint at a site if all adjoining \
links are present. *)
              If[AllTrue[Table[
                  {{dir, x}, {dir, x - id[dir]}}, {dir, 3}], 
                KeyExistsQ[index, #] &, 2],
                 gBasis += 1;
                 coords[gBasis] = x; 
                 index[x] = gBasis]],
        {i, -Floor[radius], 
         Ceiling[radius]}, {j, -Floor[radius], Ceiling[radius]},
        {k, -Floor[radius], Ceiling[radius]}]])

makeMatrixW[] := (aa = SparseArray[Flatten[{
      (* diagonal elements *) 
      Table[{i, i} -> 4, {i, basis}],
      (* iterate over plaquettes *)
      Table[Block[{x = {i, j, k}, dir1 = dirs[[1]], dir2 = dirs[[2]], 
         ii, sgn = {1, 1, -1, -1}},
        (* the four sides of the plaquette *)
        
        ii = {{dir1, x}, {dir2, x + id[dir1]}, {dir1, 
           x + id[dir2]}, {dir2, x}}; 
        Table[If[
          KeyExistsQ[index, ii[[k1]]] && KeyExistsQ[index, ii[[k2]]] &&
          k1 != k2,
          {index[ii[[k1]]], index[ii[[k2]]]} -> sgn[[k1]]*sgn[[k2]],
          Nothing],
              {k1, 4}, {k2, 4}]] ,
            {dirs, {{1, 2}, {2, 3}, {3, 1}}},
            {i, -Floor[radius], Ceiling[radius]},
            {j, -Floor[radius], Ceiling[radius]},
            {k, -Floor[radius], Ceiling[radius]}],
      (* Coulomb gauge constraints *)
      (* Scale this smaller so that Maxwell's equations
        take precidence when there is a conflict. *)
      Table[Block[{x = coords[k], gc=0.001}, 
        If[! KeyExistsQ[index, {dir, x}], 
         Print["bad key 1 ", {dir, x}]]; 
        If[! KeyExistsQ[index, {dir, x - id[dir]}], 
           Print["bad key 2 ", {dir, x - id[dir]}]];
        {{k, index[{dir, x}]} -> gc,
         {k, index[{dir, x - id[dir]}]} -> -gc}],
            {dir, 3}, {k, basis + 1, gBasis}]}],
                                   {gBasis, basis}]; 
  Block[{dir1, x, y, epsilon, kk}, 
   Do[{dir1, x} = coords[k];(* radius border *)
    (* links parallel to a given link *)
    Do[If[dir2 != dir1 && !KeyExistsQ[index, {dir1, x + sign*id[dir2]}], 
      y = x + id[dir1]/2 - id[iDir]/2; 
      If[Not[Norm[y] <= radius && Norm[y + sign*id[dir2]] >= radius], 
       Print["bad parallel link ", {{dir1, x}, dir2, sign}]; 
       Print["norms ", 
        N[{Norm[y], radius, Norm[y + sign*id[dir2]]}]]]; 
      epsilon = 
       Simplify[ Sqrt[radius^2 - y.y + y[[dir2]]^2] - sign*y[[dir2]]];
       If[epsilon < 0 || epsilon > 1, 
       Print["bad parallel epsilon ", epsilon, 
        " for ", {{dir1, x}, dir2, sign}]]; 
       aa[[k, k]] = aa[[k, k]] + 1/N[epsilon] - 1],
       {sign, {-1, 1}}, {dir2, 3}];
    (* Links perpendicular to a given link *)
    Do[kk = {dir2, x + shift1*id[dir1] + shift2*id[dir2]}; 
     If[dir2 != dir1 && 
       KeyExistsQ[index, kk] && 
       !KeyExistsQ[index,
                   {dir2, x + shift1*id[dir1] - (1 + shift2)*id[dir2]}], 
      y = kk[[2]] + id[dir2]/2 - id[iDir]/2; 
      If[Not[Norm[y] <= radius && 
         Norm[y - (2*shift2 + 1)*id[dir2]] >= radius], 
       Print["bad perp link ", {{dir1, x}, dir2, shift1, shift2}]; 
       Print["norms ", 
        N[{Norm[y], radius, Norm[y - (2*shift2 + 1)*id[dir2]]}]]]; 
      epsilon = Simplify[Sqrt[
          radius^2 - y.y + y[[dir2]]^2] + (2*shift2 + 1)*y[[dir2]]]; 
      If[epsilon < 0 || epsilon > 1, 
       Print["bad perp epsilon ", epsilon, 
        " for ", {{dir1, x}, dir2, shift1, shift2}]]; 
      aa[[k, index[kk]]] += If[
          EvenQ[Mod[dir2 - dir1, 3] + shift1 + shift2], 1, -1]*
                            (1/N[epsilon] - 1)],
       {shift1, 0, 1}, {shift2, -1, 0}, {dir2, 3}], {k, basis}]];
  b = Table[0, {gBasis}]; 
  b[[index[{iDir, {0, 0, 0}}]]] = 1.0;)


(*  Magnetostatics:  3 dimensions, using lattice symmetries *)

makeCoordsV[r_, gauge_:"Coulomb"] := (iDir = 3; basis = 0; 
  radius = r; coords = Association[];
  index = Association[];
  Do[Block[{x = {i, j, k}, dx}, dx = x + id[dir]/2 - id[iDir]/2;
    (* Exclude links with midpoints on the boundary *)
    If[Norm[dx] < radius && AllTrue[dx, # >= 0 &] && 
       Not[IntegerQ[gauge] && dir == gauge],
       basis += 1; 
       coords[basis] = {dir, x}; index[{dir, x}] = basis]],
     {dir, 3}, {i, 0, Ceiling[radius]}, {j, 0, Ceiling[radius]},
     {k, 0, Ceiling[radius]}];
  gBasis = basis; 
  If[gauge === "Coulomb",
     Do[Block[{x = {i, j, k}},
         (* Apply Coulomb gauge constraint at a site if all adjoining
           links are present in the basis.
           For sites on a reflection plane, the gauge is fixed by
           demanding the appropriate reflection symmetry. *)     
         If[AllTrue[Table[{{dir, x}, {dir, x - id[dir]}}, {dir, 3}], 
                    KeyExistsQ[index, #] &, 2],
            gBasis += 1; coords[gBasis] = x; 
            index[x] = gBasis]],
        {i, 0, Ceiling[radius]},
        {j, 0, Ceiling[radius]},
        {k, 0, Ceiling[radius]}]])

makeMatrixV[] := (aa = SparseArray[Flatten[{
      (* diagonal elements *) 
      Table[{i, i} -> 4, {i, basis}],
      (* iterate over plaquettes *)
      Table[Block[{x = {i, j, k}, dir1 = dirs[[1]], dir2 = dirs[[2]], 
         ii, sgn = {1, 1, -1, -1}},
        (* the four sides of the plaquette *)
        ii = {{dir1, x}, {dir2, x + id[dir1]},
              {dir1, x + id[dir2]}, {dir2, x}}; 
        Table[If[
          KeyExistsQ[index, ii[[k1]]] && KeyExistsQ[index, ii[[k2]]] &&
          k1 != k2,
          {index[ii[[k1]]], index[ii[[k2]]]} -> sgn[[k1]]*sgn[[k2]],
          Nothing],
              {k1, 4}, {k2, 4}]] ,
            {dirs, {{1, 2}, {2, 3}, {3, 1}}},
            {i, 0, Ceiling[radius]},
            {j, 0, Ceiling[radius]},
            {k, 0, Ceiling[radius]}],
      (* Coulomb gauge constraints *)
      (* Scale size smaller so that Maxwell's equations
        take precidence when there is a conflict. *)
      Table[Block[{x = coords[k], gc = 0.001}, 
        If[!KeyExistsQ[index, {dir, x}], 
         Print["bad key 1 ", {dir, x}]]; 
        If[!KeyExistsQ[index, {dir, x - id[dir]}], 
           Print["bad key 2 ", {dir, x - id[dir]}]];
        {{k, index[{dir, x}]} -> gc,
         {k, index[{dir, x - id[dir]}]} -> -gc}],
            {dir, 3}, {k, basis + 1, gBasis}]}],
                                   {gBasis, basis}];
 Block[{dir1, x, y, epsilon, kk}, 
  Do[{dir1, x} = coords[k];
    (* radius border *)
    (* links parallel to a given link *)
     Do[If[dir2 != dir1 && !KeyExistsQ[index, {dir1, x + id[dir2]}],
           If[dir1 == iDir,
      y = x + id[dir1]/2 - id[iDir]/2; 
      If[Not[Norm[y] <= radius && Norm[y + id[dir2]] >= radius], 
       Print["bad parallel link ", {{dir1, x}, dir2}]; 
       Print["norms ", 
        N[{Norm[y], radius, Norm[y + id[dir2]]}]]]; 
      epsilon = 
       Simplify[Sqrt[radius^2 - y.y + y[[dir2]]^2] - y[[dir2]]];
       If[epsilon < 0 || epsilon > 1, 
       Print["bad parallel epsilon ", epsilon, 
        " for ", {{dir1, x}, dir2}]]; 
       aa[[k, k]] += 1/N[epsilon] - 1,
              If[KeyExistsQ[index, {dir1, x - id[dir2]}],
              aa[[k, k]] -= 1;
              aa[[k, index[{dir1, x - id[dir2]}]]] += 1]]],
       {dir2, 3}];
    (* Links perpendicular to a given link *)
    Do[kk = {dir2, x + shift1*id[dir1] - id[dir2]}; 
     If[dir2 != dir1 && KeyExistsQ[index, kk] &&
        !KeyExistsQ[index, {dir2, x + shift1*id[dir1]}], 
      y = kk[[2]] + id[dir2]/2 - id[iDir]/2; 
      If[Not[Norm[y] <= radius && Norm[y + id[dir2]] >= radius], 
       Print["bad perp link ", {{dir1, x}, dir2, shift1}]; 
       Print["norms ", 
        N[{Norm[y], radius, Norm[y + id[dir2]]}]]]; 
      epsilon = Simplify[Sqrt[radius^2 - y.y + y[[dir2]]^2] - y[[dir2]]]; 
      If[epsilon < 0 || epsilon > 1, 
       Print["bad perp epsilon ", epsilon, 
        " for ", {{dir1, x}, dir2, shift1}]]; 
      aa[[k, index[kk]]] += If[EvenQ[Mod[dir2 - dir1, 3] + shift1], 
                               1, -1]*(1/N[epsilon] - 1)],
       {shift1, 0, 1}, {dir2, 3}];
    (* link on boundary coplanar to current, 
    assume fields are even under reflections. *)
    Do[If[dir2 != iDir && dir2 != dir1 && x[[dir2]] == 0, 
      addTo[k, {dir1, x + id[dir2]}, -1]], {dir2, 3}];
    (* link transverse to current, 
    assume fields are odd under reflections. *)
    If[iDir != dir1 && x[[iDir]] == 1, addTo[k, {dir1, x}, 1]];
    (* link on boundary traverse to current, 
    assume fields are odd under reflections. *)
    Do[If[dir1 == iDir && dir2 != dir1 && x[[iDir]] == 0, 
      addTo[k, {dir2, x + id[iDir]}, 1]; 
      addTo[k, {dir2, x - id[dir2] + id[iDir]}, -1]], {dir2, 3}],
      {k, basis}]];
  b = Table[0, {gBasis}]; 
  b[[index[{iDir, {0, 0, 0}}]]] = 1.0;)
