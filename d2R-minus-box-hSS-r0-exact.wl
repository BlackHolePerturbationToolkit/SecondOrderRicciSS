(* ::Package:: *)

(* Cancellation-free closed form of d2R[hS,hS] - Box[hSS] at r = r0, for all r0 and all alpha.

   At r = r0 the components of d2R and Box[hSS] diverge like rho^-4 at the pole (alpha = 0) while their
   difference stays finite, so evaluating them separately needs ever more precision as alpha -> 0.
   With Delta r = 0, rho = q y s exactly (s = Sin[alpha/2], q = r0 Sqrt[(r0-2)/(r0-3)], y = Sqrt[chi],
   chi = 1 - Sin[beta]^2/(r0-2)), and every alpha dependence is a rational function of s and
   c = Cos[alpha/2] plus a single Log[2 q y s]. Multiplying by (1 - s^2)^K (K the largest power of 1/c)
   turns each component into a Laurent polynomial in s (times c^0/c^1 and L^0/L^1, L = Log[s]); the
   coefficients of the negative powers of s cancel identically between d2R and Box[hSS], which leaves
       F = Sum_{k >= 0} s^k (A_k + c B_k + L (C_k + c D_k)) / (1 - s^2)^K,
   free of cancellation, valid for all 0 < alpha < Pi.

   The surds (Omega0 = w = r0^(-3/2), q, y) are symbols during the algebra and restored at the end.
   Each top-level term is expanded separately (in parallel); coefficients are summed, then reduced
   (q^2, y^2, w^2 -> rational) and simplified.

   Usage: wolframscript -file d2R-minus-box-hSS-r0-exact.wl [components] [maxTerms] [nkernels]
     components: e.g. "{1}" (default all 10); maxTerms: only the first n terms (testing; no check)
   Output: d2R-minus-box-hSS-r0-exact.m: {"K" -> {K per component}, "components" -> {tt, tr, rr, ta, tb,
   ra, rb, aa, ab, bb}, "numerators" -> {Sum_k s^k (...) per component}} in R0, b, s, c, L; F = num/(1-s^2)^K. *)

$here = DirectoryName[$InputFileName];
args = Rest[$ScriptCommandLine];
comps = If[Length[args] >= 1, ToExpression[args[[1]]], Range[10]];
maxTerms = If[Length[args] >= 2, ToExpression[args[[2]]], Infinity];
Get[FileNameJoin[{$here, "d2R.wl"}]];
body[cf_] := With[{fn = FirstCase[List @@ cf, _Function]}, fn[a, b]];
t0 = AbsoluteTime[];
d = Quiet[body /@ d2R\[Mu]\[Nu][R0, R0][[2]]] /. If[_, _, e_] :> e;
x = Quiet[body /@ BoxhSS\[Mu]\[Nu][R0, R0][[2]]] /. If[_, _, e_] :> e;
rho2 = ((-2 + R0)*R0^2*Sin[a/2]^2*(1 - Sin[b]^2/(-2 + R0)))/(-3 + R0);
desurd = {Power[R0^-3, p_Rational] :> w^(2 p), Power[rho2, p_Rational] :> (q y Sin[a/2])^(2 p)};
diffs = (d - x) /. desurd;
Clear[d, x];
(* alpha -> s, c, L *)
toSC[e_] := e /. {Log[z_] /; !FreeQ[z, a] :> (Log[z] /. Log[2 q y Sin[a/2]] -> Log[2] + lq + ly + L)} /.
    f_[arg_] /; MemberQ[{Sin, Cos, Tan, Cot, Sec, Csc}, f] && !FreeQ[arg, a] :> TrigExpand[f[arg /. a -> 2 tt]] /.
    {Sin[tt] -> s, Cos[tt] -> c, Tan[tt] -> s/c, Cot[tt] -> c/s, Sec[tt] -> 1/c, Csc[tt] -> 1/s};
cpow[e_] := Max[0, -Min[0, Cases[e, Power[c, n_Integer] :> n, Infinity], If[FreeQ[e, c], 0, -1]]];

LaunchKernels[If[Length[args] >= 3, ToExpression[args[[3]]], $ProcessorCount]];
(* one term -> Association {power of s, power of c (0/1), power of L (0/1)} -> coefficient *)
termCoefficients[t_, K_] := Module[{e = Expand[toSC[t] (1 - s^2)^K]},
   e = e /. Power[c, n_Integer] :> c^Mod[n, 2] (1 - s^2)^Quotient[n, 2];
   e = Expand[e];
   Merge[Map[{Exponent[#, s, Min] /. -Infinity -> 0, Exponent[#, c], Exponent[#, L]} -> (# /. {s -> 1, c -> 1, L -> 1}) &,
       If[Head[e] === Plus, List @@ e, {e}]], Total]];
DistributeDefinitions[toSC, termCoefficients];
rel = {q2 -> R0^2 (R0 - 2)/(R0 - 3), y2 -> 1 - Sin[b]^2/(R0 - 2), w2 -> R0^-3};
reduce[e_] := Expand[e] /. {q^n_Integer :> q^Mod[n, 2] q2^Quotient[n, 2], y^n_Integer :> y^Mod[n, 2] y2^Quotient[n, 2],
     w^n_Integer :> w^Mod[n, 2] w2^Quotient[n, 2]} /. rel;
restore = {q -> R0 Sqrt[(R0 - 2)/(R0 - 3)], y -> Sqrt[1 - Sin[b]^2/(R0 - 2)], w -> R0^(-3/2),
   lq -> Log[R0 Sqrt[(R0 - 2)/(R0 - 3)]], ly -> Log[Sqrt[1 - Sin[b]^2/(R0 - 2)]]};

result = Association[Table[
   Module[{terms, K, parts, coef, neg, t1 = AbsoluteTime[]},
     terms = If[Head[diffs[[k]]] === Plus, List @@ diffs[[k]], {diffs[[k]]}];
     If[maxTerms < Length[terms], terms = Take[terms, maxTerms]];
     K = Max[cpow[toSC[#]] & /@ terms];
     parts = With[{KK = K}, ParallelMap[termCoefficients[#, KK] &, terms, Method -> "CoarsestGrained"]];
     coef = Merge[parts, Total];
     Print["component ", k, ": ", Length[terms], " terms, K = ", K, ", ", Length[coef], " (s, c, L) powers, s from ",
       Min[Keys[coef][[All, 1]]], " to ", Max[Keys[coef][[All, 1]]], "; expanded in ", Round[AbsoluteTime[] - t1], " s"];
     coef = ParallelMap[Together[reduce[#]] &, coef];
     neg = Select[Keys[coef], #[[1]] < 0 &];
     Print["   negative powers of s: ", Length[neg], ", all vanish: ", AllTrue[neg, coef[#] === 0 &],
       "; reduced in ", Round[AbsoluteTime[] - t1], " s"];
     k -> {K, DeleteCases[KeyDropFrom[coef, neg], 0]}],
   {k, comps}]];

(* closed form, and check against the 80-digit closed forms at r0 = 17.9, beta = 6/5 *)
F[k_] := With[{K = result[k][[1]], co = result[k][[2]]},
   Total[KeyValueMap[(#2 /. restore) s^#1[[1]] c^#1[[2]] L^#1[[3]] &, co]]/(1 - s^2)^K];
If[maxTerms === Infinity,
  Do[With[{fn = FirstCase[List @@ d2R\[Mu]\[Nu][179/10, 179/10][[2, k]], _Function], fb = FirstCase[List @@ BoxhSS\[Mu]\[Nu][179/10, 179/10][[2, k]], _Function], Fk = F[k]},
     Do[Module[{cf, ex},
        cf = N[(fn[SetPrecision[al, 100], SetPrecision[6/5, 100]] - fb[SetPrecision[al, 100], SetPrecision[6/5, 100]]) /. If[_, _, e_] :> e, 30];
        ex = N[Fk /. {R0 -> 179/10, b -> 6/5, s -> Sin[al/2], c -> Cos[al/2], L -> Log[Sin[al/2]]}, 30];
        Print["  component ", k, " alpha ", N[al], ": closed forms ", N[cf, 16], "  exact form ", N[ex, 16], "  rel. diff ", ScientificForm[N[Abs[ex/cf - 1]], 2]]],
       {al, {10^-6, 10^-3, 1/10, 1, 3}}]], {k, comps}];
  Put[{"components" -> comps, "K" -> (result[#][[1]] & /@ comps), "numerators" -> ((F[#] (1 - s^2)^result[#][[1]]) & /@ comps)},
    FileNameJoin[{$here, "d2R-minus-box-hSS-r0-exact.m"}]];
  Print["wrote d2R-minus-box-hSS-r0-exact.m in ", Round[AbsoluteTime[] - t0], " s"]];
