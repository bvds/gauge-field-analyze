addChi2::usage = "Calculate the chi^2 statistic of a model.";
SetAttributes[addChi2, HoldFirst];
addChi2[model_, err2_] := 
    Block[{x = model["FitResiduals"], m = Head[model]},
          Unprotect[Evaluate[m]];
          model["chiSquared"] = If[
              MatrixQ[err2], x.LinearSolve[err2, x],
              x.(x/err2)];
          Protect[Evaluate[m]]];

covariantFit1::usage = "Wrapper for the function NonlinearModelFit that
includes a covariance matrix.  It does not handle the simple data format or \
handle constraints.  It is computationally inefficient since it evaluates \
\"form\" n^2 times for every step, where n is the number of data points.  \
Methods associated with the fit function itself will not work."; 
covariantFit1::SlotConflict = "Form contains a slot variable \"#n\" which \
may conflict with the \"Function\" method.";
Options[covariantFit1] = 
 FilterRules[
  Options[NonlinearModelFit], {Except[Weights], 
   Except[VarianceEstimatorFunction]}]; 
covariantFit1[cov_?MatrixQ, data_?MatrixQ, Except[_List, form_], pars_, vars_,
   opts : OptionsPattern[]] := 
 Block[{err2, oo, ii, transform, ff},
  If[False, Print["Dimensions: ", {Dimensions[cov], Dimensions[data]}]];
  {err2, oo} = Eigensystem[cov];
  If[Min[err2] <= 10^-10*Abs[Max[err2]],
     Message[General::npdef, cov];
     Return[$Failed]];
  transform = 
  oo.Map[(form /. MapThread[Rule, {vars, Drop[#, -1]}]) &, data];
  ff = NonlinearModelFit[oo.Map[Last, data],
    (* Expand the variable transform, 
    but don't take the part unless it is numeric. *)
    If[ii > 0, #[[ii]]]&[transform], pars, ii,
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

covariantFit2::usage = "Version of the function NonlinearModelFit \
that includes a covariance matrix for the data.  It does not handle \
the simple data format or handle constraints.  It is computationally \
efficient but only some methods of NonlinearModelFit are implemented.

The default is to use Levenberg-Marquardt for the minimization.  \
This is generally much faster, but requires calculation of the \
eigenvectors of the covariance matrix.

Alternatively, one can just solve the linear system.  \
This is generally slower, but can handle a sparse covariance matrix \
without losing sparsity.";
Options[covariantFit2] = Options[FindMinimum];
covariantFit2::singular = "Hessian is singular for `1`";
covariantFit2[cov_?MatrixQ, data_?MatrixQ, Except[_List, form_], pars_, vars_,
   opts : OptionsPattern[]] := 
 Block[{err2, oo, diff, ls, fm, fpars, jac},
  fpars = Map[If[ListQ[#], First[#], #]&, pars];
  diff = Map[(Last[#] - form /. 
             MapThread[Rule, {vars, Drop[#, -1]}]) &, data];
  fm = If[OptionValue[Method] === "LevenbergMarquardt" || 
          OptionValue[Method] === Automatic,
          {err2, oo} = Eigensystem[cov];
          If[Min[err2] <= 10^-10*Abs[Max[err2]],
             Message[General::npdef, cov];
             Return[$Failed]];
          jac = D[diff, {fpars}];
          If[False,
             Print["Jacobian ", {Dimensions[jac], pars}];
             Print["Jacobian ", MatrixForm[jac]]];
          FindMinimum[
              Null, pars, Method -> {
                  "LevenbergMarquardt",
                  (* The Sqrt[2] is from the definition of the
                    residual:  f = r.r/2.
                    The RuleDelayed allows us to evaluate the dot
                    product after the parameters have numerical
                    values. *)
                  "Residual" :> Sqrt[2/err2] * (oo.diff),
                  "Jacobian" :> Sqrt[2/err2] * (oo.jac)},
              (* FindMinimum is HoldAll *)
              Evaluate[Apply[Sequence, FilterRules[{opts}, {Method}]]]],
          ls = LinearSolve[cov];
          FindMinimum[diff.ls[diff], pars, opts]];
  (* Use Module to create a closure rather than use
    the Mathematica strategy of using subvalues; see
    https://mathematica.stackexchange.com/questions/15945 *)
  Module[{model, dof = Length[data] - Length[pars],
          (* Pre-compute values even if not used ... *)
          hess = Outer[(D[diff, #1]/.fm[[2]]).LinearSolve[
              cov, D[diff, #2]/.fm[[2]]]&, fpars, fpars]},
         Format[model] = Panel[fm];
         SetAttributes[model, Listable]; 
         model["chiSquared"] = fm[[1]];
         model["FitResiduals"] = diff/.fm[[2]];
         model["BestFitParameters"] = fm[[2]];
         model["CovarianceMatrix"] := Block[
             {ns = NullSpace[hess]},
             If[Length[ns] == 0, Inverse[hess],
                Message[covariantFit2::singular,
                        Chop[ns].Map[First, model["BestFitParameters"]]];
                $Failed]];
         model["CorrelationMatrix"] := Block[
             {covMat = model["CovarianceMatrix"], w},
             If[covMat === $Failed,
                $Failed,
                w = 1/Sqrt[Diagonal[model["CovarianceMatrix"]]];
                Map[w*#&, model["CovarianceMatrix"]]*w]];
         (* See https://mathematica.stackexchange.com/questions/89626 *)
         model["ParameterTableEntries"] := Block[
             {covMat = model["CovarianceMatrix"]},
             If[covMat === $Failed,
                Map[{#[[2]], Null, Null, Null}&,
                    model["BestFitParameters"]],
                Map[
                    Join[#, {#[[1]]/#[[2]],
                             2 (1 - CDF[StudentTDistribution[dof],
                                        Abs[#[[1]]/#[[2]]]])}]&,
                        Transpose[
                            {Map[#[[2]]&, model["BestFitParameters"]],
                             model["ParameterErrors"]}]]]];
         model["ParameterErrors"] := Block[
             {covMat = model["CovarianceMatrix"]},
             If[covMat === $Failed,
                Map[Null&, model["BestFitParameters"]],
                Sqrt[Diagonal[covMat]]]];
         model["ParameterTable"] := Style[TableForm[
             model["ParameterTableEntries"],
             TableHeadings->{Map[First, model["BestFitParameters"]],
                             {"Estimate", "Standard Error",
                              "t-Statistic", "P-Value"}}], "DialogStyle"];
         model["Function"] = Function[Evaluate[
             form/.Join[model["BestFitParameters"],
                        MapIndexed[Rule[#1, Apply[Slot, {#2[[1]]}]]&, vars]]]];
         model]];
covariantFit2[err2_?VectorQ, args__, opts:OptionsPattern[]]:=
    Block[{ff = NonlinearModelFit[args,  VarianceEstimatorFunction -> (1 &),
                                  Weights -> 1/err2, opts]},
          addChi2[ff, err2]; ff];
