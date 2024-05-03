(* ::Package:: *)

LaunchKernels[44];
Get["../d2R.wl"];
lmax = 2;
mpmax = 10;
data[r0] = ParallelTable[
  int = d2R\[Mu]\[Nu][r, r0];
  boxint = BoxhSS\[Mu]\[Nu][r,r0];
  Table[If[l==0, Print[{r0, r, i}]];
  If[r == r0, Indeterminate, d2R[int, i, l, mp]]
  , {i, 10}, {l, 0, lmax}, {mp, -Min[l,mpmax], Min[l,mpmax]}], {r, r0 - 2, r0 + 2, 1/100}, Method -> "FinestGrained"];
Save["d2R-minus-box-hSS-data-r0_" <> ToString[r0] <> ".m", data];
