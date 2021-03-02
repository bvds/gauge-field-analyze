latticeSpacing::usage = "The lattice spacing in units of the string tension for a specific beta.  From Athenodorou-Teper2017_Article_SUNGaugeTheoriesIn21Dimensions.pdf";
latticeSpacing[nd_, nc_, "1", beta_, plaquette_]:=
    latticeSpacing[nd, nc, 1, beta, plaquette];
latticeSpacing[nd_, nc_, beta_, plaquette_]:=
    latticeSpacing[nd, nc, 1, beta, plaquette];

Module[{
    (* The c1 coefficients in Table 19 in Athenodorou Teper2017 paper
      disagree with the slope of Fig. 5 and with the SU(4) "1test"
      calculation below.  Thus, we use
      c1 coefficients from arXiv:hep-lat/0206027v1 24 Jun 2002, Table 5
      Older calculation:  arXiv:hep-lat/0107007v2 1 Aug 2001, Table 11
     *)
    data = {Null,
            {valueError[0.16745, 0.00011], valueError[-0.036, 0.002],
             valueError[-0.417, 0.024]/(2 2^2)^2},
            {valueError[0.18389, 0.00017], valueError[-0.042, 0.003],
             valueError[-2.43, 0.11]/(2 3^2)^2},
            {valueError[0.18957, 0.00021], valueError[-0.038, 0.004],
             valueError[-7.74, 0.73]/(2 4^2)^2},
            {valueError[9.657, 0.038]/(2 5^2), Null,
             valueError[-21.4, 1.9]/(2 5^2)^2},
            {valueError[0.19329, 0.00011], valueError[-0.035, 0.004],
             Null},
            Null,
            {valueError[0.19486, 0.00017], valueError[-0.033, 0.005],
             Null}}},
  latticeSpacingData[3] = data;
  latticeSpacing[3, nc_, 1, beta_, plaquette_:1] :=
    Block[{betai = beta plaquette/(2 nc^2)},
          (data[[nc, 1]] + data[[nc, 3]]/betai)/betai]];


teperMass::usage = "Glueball masses in units of the lattice spacing,
  from arXiv:hep-lat/9804008v2 23 Nov 1998, Tables 25 and 26";
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

(* From Athenodorou-Teper2017_Article_SUNGaugeTheoriesIn21Dimensions.pdf
  Table 10 and 15
  Older version:  arXiv:hep-lat/0107007v2 1 Aug 2001, Table 11 *)
Module[{
    (* beta, L, plaquette, 2A and 2S and Sqrt[sigma_1]*a *)
    data = {{28, 14,  valueError[0.8093380, 0.0000007],
             valueError[1.1627, 0.0066], valueError[.99, 0.02],
             valueError[0.25205, 0.00033]},
            {40, 22, valueError[0.86961166, 0.00000024],
             valueError[0.8101, 0.0024], valueError[.42, 0.0074],
             valueError[0.16758, 0.00020]},
            {51, 30, valueError[0.89878825, 0.00000022],
             valueError[0.6519, 0.0028], valueError[.150, 0.020],
             valueError[0.12838, 0.00017]},
            {63, 36, valueError[0.91861509, 0.00000013],
             valueError[0.4941, 0.0014], valueError[0.853, 0.010],
             valueError[0.10252, 0.00011]},
            {74, 44, valueError[0.93099416, 0.00000007],
             valueError[0.4346, 0.0010], valueError[0.7623, 0.0072],
             valueError[0.08637, 0.00010]},
            {86, 50, valueError[0.94081019, 0.00000005],
             valueError[0.3578, 0.0007], valueError[0.6306, 0.0037],
             valueError[0.07378, 0.00008]}},
    tensionFromMass = Function[{ll, ee},
         (Pi/3 + Sqrt[(Pi/3)^2 + (2 ll ee)^2])/(2 ll^2)],
    d2},
   d2 = Map[{#[[1]], #[[2]], #[[3]],
             Sqrt[tensionFromMass[#[[2]], #[[4]]]],
             Sqrt[tensionFromMass[#[[2]], #[[5]]]],
             #[[6]]}&, data];
   latticeSpacingData[3, 4] = d2;
   (* May want to eventually memoize this. *)
   latticeSpacing[3, nc:4, op_, beta_, plaquette_:1] :=
    Block[{model, beta0, betai = beta*plaquette/(2 nc^2),
           col = Which[op=="2A", 4, op=="2S", 5, op=="1test", 6]},
     model = LinearModelFit[
         Map[{#[[1]]*#[[3, 1]]/(2 nc^2), #[[col, 1]]}&, d2],
         {1/beta0, 1/beta0^2}, beta0,
         IncludeConstantBasis -> False,
         VarianceEstimatorFunction -> (1 &), 
         Weights -> Map[1/#[[col, 2]]^2&, d2]];
     Print[model["ParameterTable"]];
     Print["{chi^2, dof}:  ",
           {Total[model["FitResiduals"]^2/Map[#[[col, 2]]&, d2]^2],
            Length[d2]}];
     valueError[model[betai], Norm[
         model["BasisFunctions"]*model["ParameterErrors"]
         /. beta0 -> betai]]]/;op=="2A"||op=="2S"||op=="1test"
];

(*  From Athenodorou-Teper2017_Article_SUNGaugeTheoriesIn21Dimensions.pdf
Table 20.  According to the paper, the finite a effects are negligible. *)
Module[{
    data = {
        {valueError[1.1649, 0.0011], valueError[1.5377, 0.0038]},
        Null,
        {valueError[1.2753, 0.0022], valueError[1.4961, 0.0041],
         valueError[1.3573, 0.0033], valueError[1.6861, 0.0067]}}},

  latticeSpacing[3, nc_, op_, beta_, plaquette_:1] :=
    latticeSpacing[3, nc, 1, beta, plaquette]*
    data[[nc-3, Which[op=="2A", 1, op=="2S", 2, op=="3A", 3, op="3M", 4]]]
];



(* From https://arxiv.org/abs/hep-lat/0607015
The functional form of the extrapolation in nc is not theoretically
justified but does fit the lattice data and the asymptotic value qualitatively
agrees with the one loop perturbative calculation. *)
sommerScale::usage = "The Sommer reference scale, in units of the string tension.";
sommerScale[3, _] := sommerScale[3]
sommerScale[3] := Sqrt[valueError[1.52, 0.01*1.52]];
meyerScale::usage = "Length scale marking the transition between the
perturbative regeme and the string regeme, in units of the string tension.";
meyerScale[3, nc_] := sommerScale[3]*(valueError[0.402, 0.045] +
                                  valueError[2.80, 0.20]/nc);
