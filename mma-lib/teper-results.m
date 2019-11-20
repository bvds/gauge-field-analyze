(* From arXiv:hep-lat/0206027v1 24 Jun 2002, Table 5 *)
teperTension[3, 2, beta_, plaquette_: 1] := 
  Block[{betai = beta plaquette}, {1.3405 - 0.417/betai, 
     Norm[{0.0031, 0.024/betai}]}/betai];
teperTension[3, 3, beta_, plaquette_: 1] := 
  Block[{betai = beta plaquette}, {3.3182 - 2.43/betai, 
     Norm[{0.0061, 0.11/betai}]}/betai];
teperTension[3, 4, beta_, plaquette_: 1] := 
 Block[{betai = beta plaquette}, {6.065 - 7.74/betai, 
    Norm[{0.021, 0.73/betai}]}/betai]; 
teperTension[3, 5, beta_, plaquette_: 1] := 
 Block[{betai = beta plaquette}, {9.657 - 21.4/betai, 
    Norm[{0.038, 1.9/betai}]}/betai]; 
teperTension[3, 6, beta_, plaquette_: 1] := 
 Block[{betai = beta plaquette}, {14.000 - 46.5/betai, 
    Norm[{0.038, 2.6/betai}]}/betai];

(* We are fitting a^2\[Sigma] *)
squareValue[{value_, error_}] := {value^2, 2 value error};
