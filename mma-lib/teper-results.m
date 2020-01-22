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

(* From arXiv:hep-lat/0107007v2 1 Aug 2001, Table 11 *)
Module[{
    beta,
    nd = 3,
    data = {{28, 12, 0.7152, 0.0039, 0.9718, 0.0095},
            {28, 16, 0.9937, 0.0058, 1.345, 0.014},
            {33, 16, 0.6622, 0.0029, 0.914, 0.006},
            {45, 24, 0.4974, 0.0022, 0.6815, 0.0040},
            {60, 32, 0.3571, 0.0014, 0.4888, 0.0015}},
    d2, tmodel,
    rerror = (#[[5]]/#[[3]]*Norm[{
        #[[6]]/#[[5]], #[[4]]/#[[3]]}])&},
   d2 = Map[Block[{z = Pi (nd-2)/(6*#[[2]])},
                  {#[[1]], #[[2]], #[[3]] + z, #[[4]],
                   #[[5]] + z, #[[6]]}]&, data];
   teperKstringData[3, 4] = d2;
   tmodel = LinearModelFit[
       Map[{#[[1]], #[[5]]/#[[2]]}&, d2],
       {1/beta, 1/beta^2}, beta,
       IncludeConstantBasis->False,
      VarianceEstimatorFunction -> (1 &), 
       Weights -> Map[(1/(#[[6]]/#[[2]])^2)&, d2]];
   teperTension[3, 4, beta0_, "2A"] :=
    (Print["The fit isn't very good.  Probably better to use the fundamental string tension and kepterKstringRatio[].",
           tmodel["ParameterTable"]];{tmodel[beta0], Norm[
            tmodel["BasisFunctions"]*tmodel["ParameterErrors"]
            /. beta -> beta0]});
    teperKstringRatio[3, 4, "2A"] = LinearModelFit[
        Map[{#[[1]], #[[5]]/#[[3]]} &, d2],
        {1, 1/beta0^2}, {beta0},
        VarianceEstimatorFunction -> (1 &), 
        Weights -> Map[(1/(rerror[#])^2)&, d2]]
];
