Clear[valueError];
valueError::usage = "valueError[a, b] represents value a with standard error b.";
valueError/:valueError[a_, ae_] + valueError[b_, be_] :=
     valueError[a+b, Norm[{ae, be}]];
valueError/:valueError[a_, ae_] * valueError[b_, be_] :=
    valueError[a*b, Norm[{ae*b, a*be}]];
valueError/:b_?NumericQ*valueError[a_, ae_] :=
     valueError[a*b, ae*Abs[b]];
valueError/:b_?NumericQ + valueError[a_, ae_] :=
     valueError[a+b, ae];
valueError/:valueError[a_, ae_]^b_?NumericQ :=
    valueError[a^b, ae Abs[b*a^(b-1)]]/;Abs[a]>ae;
valueError/:Exp[valueError[a_, ae_]] :=
     valueError[Exp[a], ae Exp[a]];
valueError/:Conjugate[valueError[a_, ae_]] :=
     valueError[Conjugate[a], Conjugate[ae]];

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
Format[valueError[a_?NumericQ, b_?NumericQ]] := printNumberError[a, b];
