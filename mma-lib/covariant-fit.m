addChi2::usage = "Calculate chi^2 statistic to a model.";
SetAttributes[addChi2, HoldFirst];
addChi2[model_, err2_] := 
    Block[{x = model["FitResiduals"], m = {Head[model]}},
          Apply[Unprotect, m];
          model["chiSquared"] = If[
              MatrixQ[err2], x.LinearSolve[err2, x],
              x.(x/err2)];
          Apply[Protect, m]];

covariantFit1::usage = "Wrapper for NonlinearModelFit that includes a \
covariance matrix.  It does not handle simple data format or \
handle constraints.  It is computationally inefficient since it \
evaluates the fit function n^2 times for every step, where n is the \
number of data points."; 
Options[covariantFit1] = 
 FilterRules[
  Options[NonlinearModelFit], {Except[Weights], 
   Except[VarianceEstimatorFunction]}]; 
covariantFit1[cov_, data_?MatrixQ, Except[_List, form_], pars_, vars_,
   opts : OptionsPattern[]] := 
 Block[{err2, oo, ii, transform, ff}, {err2, oo} = Eigensystem[cov]; 
  transform = 
   oo.Map[(form /. MapThread[Rule, {vars, Drop[#, -1]}]) &, data];
  ff = NonlinearModelFit[oo.Map[Last, data] ,
    (* Expand the variable transform, 
    but don't take the part unless it is numeric. *)
    If[ii > 0, #[[ii]]] &[transform], pars, ii, 
       VarianceEstimatorFunction -> (1 &), Weights -> 1/err2, opts];
  addChi2[ff, err2];
  ff];

covariantFit2::usage = "Version of NonlinearModelFit that includes a \
covariance matrix.  It does not handle the simple data format or \
handle constraints.  It is computationally efficient but only some of \
the methods of NonlinearModelFit are implemented.";
Options[covariantFit2] = Options[FindMinimum]; 
covariantFit2[cov_, data_?MatrixQ, Except[_List, form_], pars_, vars_,
   opts : OptionsPattern[]] := 
 Block[{err2, oo, diff, ls}, 
  diff = Map[(Last[#] - form /. 
       MapThread[Rule, {vars, Drop[#, -1]}]) &, data]; 
  If[OptionValue[Method] === "LevenbergMarquardt" || 
    OptionValue[Method] === Automatic, {err2, oo} = Eigensystem[cov];
   FindMinimum[Total[(oo.diff)^2/err2], pars, opts], 
   ls = LinearSolve[cov]; FindMinimum[diff.ls[diff], pars, opts]]];
