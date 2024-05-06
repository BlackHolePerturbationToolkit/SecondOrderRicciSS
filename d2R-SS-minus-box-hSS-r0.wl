(* ::Package:: *)

LaunchKernels[44];
Get["../d2R.wl"];
lmax = 2;
mpmax = 10;
data[r0] = ParallelTable[
  Table[If[l==0, Print[{r0, r, i}]];
  If[r == r0,
      int = d2R\[Mu]\[Nu][r+10^-4, r0];
      boxint = BoxhSS\[Mu]\[Nu][r+10^-4,r0];
      int2 = d2R\[Mu]\[Nu][r-10^-4, r0];
      boxint2 = BoxhSS\[Mu]\[Nu][r-10^-4,r0];
      (d2Rminusbox[int, boxint, i, l, mp]+d2Rminusbox[int2, boxint2, i, l, mp])/2
      ,
      int = d2R\[Mu]\[Nu][r, r0];
      boxint = BoxhSS\[Mu]\[Nu][r,r0];
      d2Rminusbox[int, boxint, i, l, mp]
  ];
  , {i, 10}, {l, 0, lmax}, {mp, -Min[l,mpmax], Min[l,mpmax]}], {r, r0 - 2, r0 + 2, 1/100}, Method -> "FinestGrained"];
Save["d2R-minus-box-hSS-data-r0_" <> ToString[NumberForm[r0,{Infinity,2}]] <> ".m", data];
