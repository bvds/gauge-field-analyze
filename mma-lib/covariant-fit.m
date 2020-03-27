addChi2::usage = "Calculate chi^2 statistic to a model.";
SetAttributes[addChi2, HoldFirst];
addChi2[model_, err2_] := 
    Block[{x = model["FitResiduals"], m = Head[model]},
          Unprotect[Evaluate[m]];
          model["chiSquared"] = If[
              MatrixQ[err2], x.LinearSolve[err2, x],
              x.(x/err2)];
          Protect[Evaluate[m]]];

covariantFit1::usage = "Wrapper for NonlinearModelFit that includes a \
covariance matrix.  It does not handle the simple data format or \
handle constraints.  It is computationally inefficient since it \
evaluates \"form\" n^2 times for every step, where n is the \
number of data points.  Methods associated with the fit function itself \
will not work."; 
covariantFit1::SlotConflict = "Form contains a slot variable \"#n\" which may conflict with the \"Function\" method.";
Options[covariantFit1] = 
 FilterRules[
  Options[NonlinearModelFit], {Except[Weights], 
   Except[VarianceEstimatorFunction]}]; 
covariantFit1[cov_?MatrixQ, data_?MatrixQ, Except[_List, form_], pars_, vars_,
   opts : OptionsPattern[]] := 
 Block[{err2, oo, ii, transform, ff},
  If[False, Print["Dimensions: ", {Dimensions[cov], Dimensions[data]}]];
  {err2, oo} = Eigensystem[cov]; 
  transform = 
   oo.Map[(form /. MapThread[Rule, {vars, Drop[#, -1]}]) &, data];
  ff = NonlinearModelFit[oo.Map[Last, data] ,
    (* Expand the variable transform, 
    but don't take the part unless it is numeric. *)
    If[ii > 0, #[[ii]]] &[transform], pars, ii, 
       VarianceEstimatorFunction -> (1 &), Weights -> 1/err2, opts];
  addChi2[ff, err2];
  Unprotect[FittedModel];
  ff["BestFit"] = form/.ff["BestFitParameters"];
  If[Not[FreeQ[form, Slot[_]]],
     Message[covariantFit1::SlotConflict]];
  ff["Function"] = Function[Evaluate[
      form/.Join[ff["BestFitParameters"],
                 MapIndexed[Rule[#1, Apply[Slot, {#2[[1]]}]]&, vars]]]];
  Protect[FittedModel];
  ff];
covariantFit1[err2_?VectorQ, args__, opts:OptionsPattern[]]:=
    Block[{ff = NonlinearModelFit[args,  VarianceEstimatorFunction -> (1 &),
                                  Weights -> 1/err2, opts]},
          addChi2[ff, err2]; ff];

covariantFit2::usage = "Version of NonlinearModelFit that includes a \
covariance matrix.  It does not handle the simple data format or \
handle constraints.  It is computationally efficient but only some of \
the methods of NonlinearModelFit are implemented.";
Options[covariantFit2] = Options[FindMinimum]; 
covariantFit2[cov_?MatrixQ, data_?MatrixQ, Except[_List, form_], pars_, vars_,
   opts : OptionsPattern[]] := 
 Block[{err2, oo, diff, ls, fm}, 
  diff = Map[(Last[#] - form /. 
       MapThread[Rule, {vars, Drop[#, -1]}]) &, data]; 
  fm = If[OptionValue[Method] === "LevenbergMarquardt" || 
    OptionValue[Method] === Automatic, {err2, oo} = Eigensystem[cov];
   FindMinimum[Total[(oo.diff)^2/err2], pars, opts], 
     ls = LinearSolve[cov];
     FindMinimum[diff.ls[diff], pars, opts]];
  (* Use Module to create a closure rather than
    the Mathematica strategy of using subvalues; see
    https://mathematica.stackexchange.com/questions/15945 *)
  Module[{model},
         Format[model] = Panel[fm];
         model["chiSquared"] = fm[[1]];
         model["BestFitParameters"] = fm[[2]];
         model[args__?NumericQ] := form/.Join[
             fm[[2]], MapThread[Rule, {vars, {args}}]];
         model]];
covariantFit2[err2_?VectorQ, args__, opts:OptionsPattern[]]:=
    Block[{ff = NonlinearModelFit[args,  VarianceEstimatorFunction -> (1 &),
                                  Weights -> 1/err2, opts]},
          addChi2[ff, err2]; ff];
