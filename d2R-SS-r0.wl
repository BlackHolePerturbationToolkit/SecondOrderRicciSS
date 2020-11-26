LaunchKernels[44];
Get["../d2R.wl"];
mpmax = 10;
data[r0] = ParallelTable[
  int = d2R\[Mu]\[Nu][r, r0];
  Table[If[l==0, Print[{r0, r, i}]];
  If[r == r0, Indeterminate, d2R[int, i, l, mp]]
  , {i, 10}, {l, 0, 20}, {mp, -Min[l,mpmax], Min[l,mpmax]}], {r, r0 - 2, r0 + 2, 1/100}, Method -> "FinestGrained"];
Save["d2R-data-r0_" <> ToString[r0] <> ".m", data];
