condition::usage = "Calculate 2-norm condition number of a matrix \
using the SVD."; 
condition[m_] := 
 If[Length[m] < 1000, (First[#]/Last[#]) &[
   SingularValueList[Normal[m], Tolerance -> 0]], 
  First[SingularValueList[m, 1, Method -> "Arnoldi"]/ 
	SingularValueList[m, -1, Method -> {"Arnoldi"}]]];
