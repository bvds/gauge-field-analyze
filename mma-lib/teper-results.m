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

(* Glueball masses in units of the lattice spacing,
  from arXiv:hep-lat/9804008v2 23 Nov 1998, Tables 25 and 26 *)
teperMass[3, 3, 28, "0++"] = 0.5517; 
teperMass[3, 3, beta0_, "0++"] = Module[{beta}, 
  Fit[{{15, 1.095}, {21, 0.7561}, {28, 0.5517}, {34, 0.4482}}, {1, 
     1/beta^2}, beta] /. beta -> beta0]; 
teperMass[3, 3, beta0_, "0--"] = Module[{beta}, 
  Fit[{{15, 1.569}, {21, 1.101}, {28, 0.8133}, {34, 0.6682}}, {1, 
     1/beta^2}, beta] /. beta -> beta0]; 
teperMass[3, 4, beta0_, "0++"] = Module[{beta}, 
  Fit[{{28, 1.083}, {40, 0.7109}, {51, 0.5466}}, {1, 1/beta^2}, 
    beta] /. beta -> beta0]; 
teperMass[3, 4, beta0_, "0--"] = Module[{beta}, 
  Fit[{{28, 1.599}, {40, 1.039}, {51, 0.8040}}, {1, 1/beta^2}, 
    beta] /. beta -> beta0];
