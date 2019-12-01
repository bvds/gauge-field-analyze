(* Do a link at a time.  This is mostly for posterity & debugging
since the resulting action density was dominated by one-plaquette spikes. *)

linkSaddlePointStep::usage = "One iteration of the saddle point \
search.  Set rescaleCutoff to a large number to test this against the \
explicit calculation."; 
Options[linkSaddlePointStep] = {rescaleCutoff -> 1, 
  dampingFactor -> 1, printHessian -> False, printShift -> False, 
  makeHess0 -> False}; 
linkSaddlePointStep[u2Staples_, u1_, OptionsPattern[]] := 
 Block[{zzz = OptionValue[rescaleCutoff], 
   cutoff = Min[0.1, If[NumberQ[beta], Sqrt[nc^2/(3 beta)], 1]], 
   damping = OptionValue[dampingFactor], debug = False, 
   ff = u2Staples.u1, grad, mm, oo, values, deltas, sqrtU, 
   nc = Length[u2Staples]}, If[False, Print["cutoff:  ", cutoff]]; 
  grad = Map[Tr[ff.#]&, SUGenerators[]];
  If[debug, Print["grad:  ", grad]];
  (* Solve linear system mm.shifts = -Im[grad], 
  paying careful attention to nullspace of mm.  
  Could equivalently use LDL^T factorization with symmetric pivoting.  *)
  mm = Re[Tr[ff]]/(2.0 Length[ff]) IdentityMatrix[Length[grad]] + 
       SUSymmetric[].Re[grad]/2.0;
  {values, oo} = Eigensystem[mm];
  (* Using the fact that the periodicity of B_a is 4 Pi or less. *) 
  deltas = -damping applyCutoff1[values, oo.Im[grad], cutoff, zzz];
  If[OptionValue[printHessian],
     Print["Hessian:  ", mm]; 
     Print["gradient:  ", Im[grad]]]; 
  If[OptionValue[printShift],
     Print["shifts:  ", deltas.oo]]; 
  If[OptionValue[makeHess0],
     gradient0 = Join[gradient0, -Im[grad]]; 
     hess0 = If[MatrixQ[hess0], 
		ArrayFlatten[{{hess0, 0}, {0, -mm}}], -mm]; 
     delta0 = Join[delta0, deltas.oo]];
  {u1.MatrixExp[I deltas.oo.SUGenerators[]], Max[Abs[deltas]]}]; 

linkSaddlePoint::usage = "Saddle point search for a single link.";
Options[linkSaddlePoint] = 
  Join[{tolerance -> 10^-6, maxCount -> 1, printHessian -> False, 
    printShift -> False, printDebug -> False, 
    updateGaugeField -> True}, Options[linkSaddlePointStep]];
linkSaddlePoint[dir_, coords_, opts : OptionsPattern[]] := 
 Block[{maxShift = 4 Pi, u1 = SUPower[getLink[dir, coords], 0.5], u2, 
   u2SumStaples, count = 0, uu}, u2 = u1;
  u2SumStaples = u2.sumStaples[dir, coords];
  While[maxShift > OptionValue[tolerance] && 
	count < OptionValue[maxCount],
	count++;
	{u1, maxShift} = linkSaddlePointStep[
	    u2SumStaples, u1, Apply[Sequence, 
	     FilterRules[{opts}, Options[linkSaddlePointStep]]]]; 
	If[OptionValue[printDebug], Print["Max shift:  ", maxShift]]]; 
  uu = u1.u2;
  If[OptionValue[printDebug], SUMatrixQ[uu]];
  If[OptionValue[updateGaugeField],
     setLink[dir, coords, uu]]; {uu, maxShift}]; 
Options[latticeSaddlePoint] = 
 Join[{maxCount -> 1, tolerance -> 10^-6}, Options[linkSaddlePoint]]; 
latticeSaddlePoint[opt:OptionsPattern[]] := 
 Block[{count = 0, maxShift = 4 Pi, linkShift, stats = {}, stat, uu, 
	newGaugeField}, 
  While[count < OptionValue[maxCount] && 
	maxShift > OptionValue[tolerance],
	maxShift = 0; 
	If[OptionValue[makeHess0],
	   gradient0 = {}; hess0 = Null; 
	   delta0 = {}]; 
	newGaugeField = 
	Table[{uu, linkShift} = 
	      linkSaddlePoint[
		  dir, latticeCoordinates[k], maxCount -> 1, 
		  Apply[Sequence,
			FilterRules[{opt},
				    Options[linkSaddlePoint]]]];
	      maxShift = Max[linkShift, maxShift]; 
	      uu,
	      {dir, nd}, {k, latticeVolume[]}]; 
	Print[stat = {count, maxShift, averagePlaquette[]}]; 
	If[False, 
	   AppendTo[stat, exponentialModel[
	       talliesToAverageErrors[polyakovLoopTallies["simple"]], 
	       printResult -> True]]];
	AppendTo[stats, stat]; 
	count++];
  {newGaugeField, stats}];
