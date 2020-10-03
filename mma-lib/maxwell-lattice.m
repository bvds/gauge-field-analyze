(* General utilities *)

id[x_] := IdentityMatrix[3][[x]];
(* uses globally defined association "index" *)
addTo[k_?IntegerQ, co_?ListQ, dx_: 1] := 
  If[KeyExistsQ[index, co], aa[[k, index[co]]] += dx];


(* Poisson equation:  2 dimensions *)

makeCoords2[r_] := (maxx = 0; basis = 0; radius = r; 
  coords = Association[];
  index = Association[]; 
  Do[If[Norm[{i, j}] < r,
        basis += 1; If[maxx < i, maxx = i]; 
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
    If[!KeyExistsQ[index, {x + 1, y}], 
     epsilon = N[Sqrt[radius^2 - y^2] - x]; 
     If[False, Print[{{x, y}, epsilon}]]; aa[[k, k]] += 1/epsilon - 1;
      If[!KeyExistsQ[index, {x, y + 1}], 
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
           If[KeyExistsQ[index, {x, y, z} + delta],
              {i, index[{x, y, z} + delta]} -> -1, 
              Nothing],
           {delta, {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
                    {0, 1, 0}, {0, 0, -1}, {0, 0, 1}}}],
              {i, i} -> 6],
                   {i, basis}]],
      basis]]; 
  Block[{x, y, z, epsilon}, 
   Do[{x, y, z} = coords[k];
    (* radius border *)
    If[!KeyExistsQ[index, {x + 1, y, z}], 
     epsilon = N[Sqrt[radius^2 - y^2 - z^2] - x]; 
     aa[[k, k]] += 1/epsilon - 1; 
     If[!KeyExistsQ[index, {x, y + 1, z}], 
      epsilon = N[Sqrt[radius^2 - x^2 - z^2] - y]; 
      aa[[k, k]] += 1/epsilon - 1]; 
     If[!KeyExistsQ[index, {x, y, z + 1}], 
      epsilon = N[Sqrt[radius^2 - x^2 - y^2] - z]; 
      aa[[k, k]] += 1/epsilon - 1]];
    (* top face *)
    If[y == z, addTo[k, {x, y, z - 1}, -1]; addTo[k, {x, y + 1, z}, -1]];
    (* bottom face *)
    If[z == 0 , addTo[k, {x, y, 1}, -1]];
    (* side face *)  
    If[x == y, addTo[k, {x + 1, y, z}, -1]; addTo[k, {x, y - 1, z}, -1]];
    (* Edge *)
    If[y == z == 0, addTo[k, {x, 1, 0}, -2]]; 
    If[z == x == y, addTo[k, {x, y, z - 1}, -1]; addTo[k, {x + 1, y, z}, -1]];
    If[x == 0 && y == 0 && z == 0, addTo[k, {1, 0, 0}, -3]],
      {k, basis}]];
  b = Table[0, {basis}];
  b[[index[{0, 0, 0}]]] = 1.0;)
