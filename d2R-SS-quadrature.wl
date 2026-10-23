(* ::Package:: *)

(* SS Ricci modes (and the matching hSS modes) by a fixed tensor-product quadrature, an alternative
   to the per-mode NIntegrate of d2R-SS-r0.wl, d2R-SS-minus-box-hSS-r0.wl and hSS-r0.wl/drhSS-r0.wl.

   The per-mode NIntegrate scatters by ~1e-5 ... 1e-4 (absolute, rotated m' = 0 modes of d2R; about
   1 % of the regular remainder S = Src_SS - E0[h2P_SS]) and occasionally fails outright (r0 = 17.9,
   r - r0 = +0.06: 4e-3): the adaptive subdivision does not always resolve the peak of width
   ~|r - r0|/r0 at the pole, and the machine-precision integrand loses digits there (roundoff
   ~1/alpha^2; the series used for alpha < 10^(-Precision/3) is truncated in r - r0 as well).

   Here, for each r:
     alpha: Gauss-Legendre panels graded geometrically from a/8 to Pi (a = |r - r0|/r0);
     beta:  trapezoid rule (periodic integrand, exponential convergence);
     the 10 tensor components are evaluated once per node and reused for every (i, l, m'):
       alpha >= alphaHP: the compiled machine-precision functions;
       alpha <  alphaHP: the uncompiled body of the same functions (exact r, r0; the series
         switch then always takes the closed form), with common subexpressions eliminated and
         constants converted to numbers, at WP digits (~10 ms per node for all 10 components).
   The mode integrands and parity rules are those of d2R[] (identical in hSS[]).

   Modes ("what"):
     "d2R"      d2R[hS,hS]                         -> src.h5 layout (Convert-d2RSS.nb), r0 rows zero
     "minusbox" d2R[hS,hS] - Box[hSS] - GD[hSS]: Box (d2R.wl/hSS.wl BoxhSS) is the plain Lorenz operator
                and GD (hSS.wl GDhSS) its gauge-damping term, together the gauge-damped Lorenz operator of
                the field equations (checked against SecondOrderSelfForce E0 with GaugeDamped -> True)
                                                    -> src.h5 layout; r0 rows evaluated at r0 itself
     "minusboxngd" d2R[hS,hS] - Box[hSS]: as minusbox without the gauge-damping term (2026-10-23), which the
                effective source can add mode by mode from the hSS modes (it is local in each mode and regular at
                r0); without it every component keeps a definite parity in beta and the other-parity modes
                (i = 2, 4, 8) vanish            -> src.h5 layout; r0 rows evaluated at r0 itself
     "gd"       GD[hSS] alone at r = r0 only (2026-10-23): the r0 rows for the mode-wise gauge-damping term of the
                minusboxngd source (lorenz2sf/ssgd.py computes it from the hSS modes elsewhere; the hSS modes are
                not available at r0) -> src.h5 layout with r = {r0, r0}
     "hSS"      hSS and d/dr hSS (hSS.wl), BLS fields, d/dr including the r dependence of the BLS
                factors (Convert-hSS.nb)            -> hP.h5 layout: "SS/i{i}l{l}m{m}" with columns
                                                       Re h, Im h, Re dh, Im dh, and "r"; r0 rows zero
   Rotation to (l, m): WignerD[{l, m, m'}, Pi, Pi/2, Pi/2], summed over all |m'| <= min(l, mpmax).

   Usage: wolframscript -file d2R-SS-quadrature.wl what r0 out.h5 [lmax] [drmax] [nkernels]
     r0 exact (e.g. 179/10); drmax: half-width of the radial range (2, or 4 for r0 > 20).
   Validated ("d2R", r0 = 17.9) against an independent refinement (30 Gauss points, 64 beta
   points, alphaHP = 1/10, 40 digits): agreement to ~1e-15 relative in the rotated modes. *)

$here = DirectoryName[$InputFileName];
(* arguments after the script name: $ScriptCommandLine under wolframscript, $CommandLine under
   WolframKernel -script *)
argv = If[Length[$ScriptCommandLine] > 1, Rest[$ScriptCommandLine],
   Module[{k = FirstPosition[$CommandLine, s_String /; StringEndsQ[s, "d2R-SS-quadrature.wl"], {0}, {1}][[1]]}, $CommandLine[[k + 1 ;;]]]];
If[Length[argv] < 3, Print["usage: d2R-SS-quadrature.wl what r0 out.h5 [lmax] [drmax] [nkernels]"]; Exit[1]];
{what, r0, out} = {argv[[1]], ToExpression[argv[[2]]], argv[[3]]};
lmax = If[Length[argv] >= 4, ToExpression[argv[[4]]], 10];
drmax = If[Length[argv] >= 5, ToExpression[argv[[5]]], 2];
nk = If[Length[argv] >= 6, ToExpression[argv[[6]]], 12];
If[!MemberQ[{"d2R", "minusbox", "minusboxngd", "gd", "hSS"}, what], Print["what must be d2R, minusbox, minusboxngd, gd or hSS"]; Exit[1]];
mpmax = 10; nGL = 20; nBeta = 32; alphaHP = 1/50; WP = 32; minDigits = 14; minAccuracy = 16; maxDoublings = 3;
(* quadrature resolution overrides (convergence checks): SS_QUAD_NGL, SS_QUAD_NBETA *)
With[{v = Environment["SS_QUAD_NGL"]}, If[StringQ[v], nGL = ToExpression[v]]];
(* geometric ratio of the alpha panels (default 2) and a test mode that evaluates only the given offsets
   r - r0 (e.g. "{0, 1/100}"), prints the time per radius and writes the rotated modes with Put *)
gradeRatio = 2; With[{v = Environment["SS_QUAD_GRADE"]}, If[StringQ[v], gradeRatio = ToExpression[v]]];
testOffsets = With[{v = Environment["SS_QUAD_TEST"]}, If[StringQ[v], ToExpression[v], None]];
With[{v = Environment["SS_QUAD_NBETA"]}, If[StringQ[v], nBeta = ToExpression[v]]];


(* ::Section:: *)
(*Definitions (loaded on every kernel)*)


defs := (
  Needs["NumericalDifferentialEquationAnalysis`"];
  (* d2R.wl: d2R and Box[hSS]; hSS.wl: hSS and the gauge-damping term GDhSS. Their small-alpha series
     branches (truncated in r - r0) are never used: the uncompiled bodies take the closed form, and the
     compiled functions are only used for alpha >= alphaHP *)
  Get[FileNameJoin[{$here, "d2R.wl"}]]; Get[FileNameJoin[{$here, "hSS.wl"}]];
  lines = StringSplit[ReadString[FileNameJoin[{$here, "d2R.wl"}]], "\n"];
  (* the mode integrands and parity rule of d2R[], with the components as symbols c[1] ... c[10] *)
  Module[{k0 = First@FirstPosition[lines, s_String /; StringStartsQ[s, "d2R[int_,i_Integer,l_Integer,mp_Integer]:="]], sw, cond, names},
    sw = StringRiffle[lines[[k0 + 7 ;; k0 + 19]], "\n"];
    cond = lines[[k0 + 22]];
    names = {"d2Rtt", "d2Rtr", "d2Rrr", "d2Rt\\[Alpha]", "d2Rt\\[Beta]", "d2Rr\\[Alpha]", "d2Rr\\[Beta]", "d2R\\[Alpha]\\[Alpha]", "d2R\\[Alpha]\\[Beta]", "d2R\\[Beta]\\[Beta]"};
    sw = StringReplace[sw, Table[names[[k]] <> "[\\[Alpha],\\[Beta]]" -> "c[" <> ToString[k] <> "]", {k, 10}]];
    sw = StringReplace[sw, {"integrand=Switch[i," -> "Switch[i,", "$Failed];" -> "$Failed]"}];
    ToExpression["modeIntegrand[i_,l_,mp_,r_,f_]:=" <> sw];
    ToExpression["modeNonzero[i_,l_,mp_]:=" <> StringTrim[cond, ","]]];
  (* the component functions (compiled, machine precision) of each quantity *)
  componentSets[r_, r0_] := Switch[what,
    "d2R", {d2R\[Mu]\[Nu][r, r0][[2]]},
    "minusbox", {d2R\[Mu]\[Nu][r, r0][[2]], BoxhSS\[Mu]\[Nu][r, r0][[2]], GDhSS\[Mu]\[Nu][r, r0][[2]]},
    "minusboxngd", {d2R\[Mu]\[Nu][r, r0][[2]], BoxhSS\[Mu]\[Nu][r, r0][[2]]},
    "gd", {GDhSS\[Mu]\[Nu][r, r0][[2]]},
    "hSS", {hSS\[Mu]\[Nu][r, r0][[2]], drhSS\[Mu]\[Nu][r, r0][[2]]}];
  combine[sets_] := Switch[what, "minusbox", {sets[[1]] - sets[[2]] - sets[[3]]}, "minusboxngd", {sets[[1]] - sets[[2]]}, _, sets];   (* {quantity(, its d/dr)} *)
  (* WP-digit evaluator of all components of all sets from the uncompiled bodies of the compiled
     functions (exact r, r0): CSE, then numeric constants (exact rationals and surds make the
     uncompiled evaluation ~400x slower) *)
  hpOpt[r_, r0_] := Experimental`OptimizeExpression[
      Map[With[{fn = FirstCase[List @@ #, _Function]}, fn[hpA, hpB] /. If[_, _, e_] :> e] &, componentSets[r, r0], {2}], "OptimizationLevel" -> 2];
  hpFunction[opt_, wp_] := Module[{held},
    held = Extract[opt, {1}, Hold] /. {q_Rational :> RuleCondition[N[q, wp + 8]],
        Power[x_Integer | x_Rational, y_Rational] :> RuleCondition[N[x^y, wp + 8]], Pi :> RuleCondition[N[Pi, wp + 8]]};
    Function @@ Join[Hold[{hpA, hpB}], held]];
  (* adaptive working precision: d2R and Box[hSS] grow like rho^-4 near the pole while d2R - Box[hSS]
     stays finite, so the digits lost to cancellation grow without bound as alpha -> 0 (all the more at
     r = r0); re-evaluate at doubled precision until every component has minDigits correct digits or
     (for components that cancel to ~0) an absolute accuracy of 10^-minAccuracy *)
  hpEval[hps_, x_, b_] := Module[{k = 1, v, wp},
    While[True,
      wp = WP 2^(k - 1);
      v = combine[hps[[k]][SetPrecision[x, wp], SetPrecision[b, wp]]];
      If[AllTrue[Flatten[v], Precision[#] >= minDigits || Accuracy[#] >= minAccuracy &] || k == Length[hps], Return[v]];
      k++]];
  (* signs of the components under beta -> -beta ({quantity, component} array of +-1), from the compiled
     functions at two generic nodes; None unless every ratio is +-1 to 1e-10 *)
  parity[cc_] := Module[{nodes = Tuples[{{0.37, 1.91}, {0.53, 2.29}}], v, w, s},
    v = Table[combine[Map[#[n[[1]], n[[2]]] &, cc, {2}]], {n, nodes}];
    w = Table[combine[Map[#[n[[1]], -n[[2]]] &, cc, {2}]], {n, nodes}];
    s = Quiet[Round[w[[1]]/v[[1]]]];
    If[!ArrayQ[s, _, MemberQ[{1, -1}, #] &], Return[None]];
    If[Max[Abs[Flatten[w - (s # &) /@ v]]] <= 10^-10 Max[Abs[Flatten[v]]], s, None]];
  glNodes[n_] := glNodes[n] = Transpose[NumericalDifferentialEquationAnalysis`GaussianQuadratureWeights[n, -1, 1, 40]];
  (* rotated modes {i, l, m', value} of each quantity at one r *)
  rotatedModes[r_, r0_] := Module[{a, edges, xs, ws, alphas, walpha, betas, wbeta, hp, cc, comps, f, A, B, W, sg, half},
    a = If[r == r0, 10^-10, Abs[r - r0]/r0];   (* r0 (minusbox, minusboxngd): no particle scale, grade from 1e-10 *)
    edges = Join[{0}, Select[a gradeRatio^Range[-3 Log[2, gradeRatio]^-1, 60 Log[2, gradeRatio]^-1], # < Pi &], {Pi}];
    {xs, ws} = glNodes[nGL];
    alphas = Flatten[Table[(edges[[k]] + edges[[k + 1]])/2 + (edges[[k + 1]] - edges[[k]])/2 xs, {k, Length[edges] - 1}]];
    walpha = Flatten[Table[(edges[[k + 1]] - edges[[k]])/2 ws, {k, Length[edges] - 1}]];
    betas = -Pi + 2 Pi Range[0, nBeta - 1]/nBeta; wbeta = 2 Pi/nBeta;
    hp = With[{opt = hpOpt[r, r0]}, Table[hpFunction[opt, WP 2^k], {k, 0, maxDoublings}]];
    cc = componentSets[N[r], N[r0]];
    (* every component is even or odd in beta (signs measured here, not assumed): evaluate at beta in [0, Pi]
       (node 1 = -Pi = Pi and nodes nBeta/2 + 1 ... nBeta) and mirror the rest; all nodes if a sign is unclear *)
    sg = parity[cc];
    half = If[sg === None, Range[nBeta], Join[{1}, Range[nBeta/2 + 1, nBeta]]];
    comps = Table[
      If[al < alphaHP,
        Table[hpEval[hp, al, b], {b, betas[[half]]}],
        With[{x = N[al]}, Table[combine[Map[#[x, N[b]] &, cc, {2}]], {b, betas[[half]]}]]],
      {al, alphas}];
    If[sg =!= None,      (* node j <-> -beta at nBeta + 2 - j *)
      comps = Map[Function[row, Table[Which[j == 1, row[[1]], j > nBeta/2, row[[j - nBeta/2 + 1]],
           True, sg row[[nBeta + 2 - j - nBeta/2 + 1]]], {j, nBeta}]], comps]];
    (* {quantity, component, alpha, beta}; machine reals (N[] would keep the WP-digit values) *)
    comps = Developer`ToPackedArray[Transpose[SetPrecision[comps, MachinePrecision], {3, 4, 1, 2}]];
    f = N[1 - 2/r];
    (* machine reals throughout the mode sums (N[] leaves the 40-digit Gauss nodes at 40 digits) *)
    A = Transpose[ConstantArray[SetPrecision[alphas, MachinePrecision], Length[betas]]];
    B = ConstantArray[SetPrecision[betas, MachinePrecision], Length[alphas]];
    W = Outer[Times, SetPrecision[walpha, MachinePrecision], ConstantArray[SetPrecision[wbeta, MachinePrecision], Length[betas]]];
    Table[{i, l, mp, If[what === "minusbox" || what === "gd" || modeNonzero[i, l, mp],   (* GD: other parity in i = 2, 4, 8 *)
        Total[W (modeIntegrand[i, l, mp, N[r], f] /. {c[k_] :> q[[k]], \[Alpha] -> A, \[Beta] -> B}), 2], 0]},
      {q, comps}, {i, 10}, {l, 0, lmax}, {mp, -Min[l, mpmax], Min[l, mpmax]}]];
  );


(* ::Section:: *)
(*Run*)


LaunchKernels[nk];
DistributeDefinitions[defs];
With[{v = {$here, what, mpmax, nGL, nBeta, alphaHP, WP, minDigits, minAccuracy, maxDoublings, lmax, gradeRatio}},
  ParallelEvaluate[{$here, what, mpmax, nGL, nBeta, alphaHP, WP, minDigits, minAccuracy, maxDoublings, lmax, gradeRatio} = v; defs, DistributedContexts -> None]];
defs;
If[testOffsets =!= None,
  Module[{res = ParallelMap[Module[{t = AbsoluteTime[], v = rotatedModes[r0 + #, r0]}, {#, AbsoluteTime[] - t, v}] &, testOffsets,
       Method -> "FinestGrained", DistributedContexts -> None]},
    Do[Print["  r - r0 = ", x[[1]], ": ", Round[x[[2]], 0.1], " s"], {x, res}];
    Put[res, out]; Print["wrote ", out]; Exit[0]]];
rs = Select[r0 + Range[-100 drmax, 100 drmax]/100, # != r0 &];
(* minusbox is regular at r0: evaluated there directly (the cluster script averaged r0 -+ 1e-4) *)
If[what === "minusbox" || what === "minusboxngd", rs = Append[rs, r0]];
(* gd: only r0 (both r0 rows of the src.h5 layout) *)
If[what === "gd", rs = {r0}];
Print[what, ", r0 = ", r0, ": ", Length[rs], " radii, lmax = ", lmax, ", ", nk, " kernels"];
t0 = AbsoluteTime[];
data = Association[ParallelMap[(# -> rotatedModes[#, r0]) &, rs, Method -> "FinestGrained", DistributedContexts -> None]];
Print["modes done in ", Round[AbsoluteTime[] - t0], " s"];

lmin[i_] := Switch[i, 1 | 2 | 3 | 6, 0, 4 | 5 | 8 | 9, 1, 7 | 10, 2];
rside = Select[r0 + Range[-100 drmax, 100 drmax]/100, # != r0 &];
rgrid = If[what === "gd", {r0, r0}, Join[Select[rside, # < r0 &], {r0, r0}, Select[rside, # > r0 &]]];
rot[q_, i_, l_, m_, rr_] := Sum[N[WignerD[{l, m, mp}, Pi, Pi/2, Pi/2]] data[rr][[q, i, l + 1, mp + Min[l, mpmax] + 1, 4]], {mp, -Min[l, mpmax], Min[l, mpmax]}];
atr0[q_, i_, l_, m_] := If[what === "minusbox" || what === "minusboxngd" || what === "gd", rot[q, i, l, m, r0], 0.];
lm[q_, i_, l_, m_] := Table[If[rr == r0, atr0[q, i, l, m], rot[q, i, l, m, rr]], {rr, rgrid}];
srcName[i_, l_, m_] := "src i=" <> ToString[i] <> " l=" <> ToString[l] <> " m=" <> ToString[m];
modesLM = Join[Flatten[Table[{i, l, m}, {i, 1, 7}, {l, lmin[i], lmax}, {m, l, 0, -2}], 2],
   Flatten[Table[{i, l, m}, {i, 8, 10}, {l, lmin[i], lmax}, {m, l - 1, 0, -2}], 2]];

Which[
  what === "d2R" || what === "minusbox" || what === "minusboxngd" || what === "gd",
  ds = Join[{"r" -> N[rgrid]}, (srcName @@ # -> ReIm[lm[1, Sequence @@ #]]) & /@ modesLM],

  what === "hSS",
  (* d/dr of the BLS fields: d/dr of the components inserted in the BLS combinations (dh) plus the
     r dependence of the BLS factors (Convert-hSS.nb, M = 1) *)
  Module[{h, dh, rv = N[rgrid]},
    Dr[1, l_, m_] := dh[1, l, m] + h[1, l, m]/(-2 + rv) - 2 h[3, l, m]/rv^2;
    Dr[2, l_, m_] := dh[2, l, m] + h[2, l, m]/(-2 + rv);
    Dr[3, l_, m_] := (rv (-2 + rv)^2 dh[3, l, m] - 2 rv h[1, l, m] + (-2 + rv)^2 h[3, l, m])/(rv (-2 + rv)^2);
    Dr[4, l_, m_] := dh[4, l, m];
    Dr[5, l_, m_] := dh[5, l, m] + 2 h[5, l, m]/((-2 + rv) rv);
    Dr[6, l_, m_] := dh[6, l, m] - h[6, l, m]/rv;
    Dr[7, l_, m_] := dh[7, l, m] - h[7, l, m]/rv;
    Dr[8, l_, m_] := dh[8, l, m];
    Dr[9, l_, m_] := dh[9, l, m] + 2 h[9, l, m]/((-2 + rv) rv);
    Dr[10, l_, m_] := dh[10, l, m] - h[10, l, m]/rv;
    h[i_, l_, m_] := h[i, l, m] = lm[1, i, l, m];
    dh[i_, l_, m_] := dh[i, l, m] = lm[2, i, l, m];
    ds = Join[{"r" -> rv}, ("SS/i" <> ToString[#[[1]]] <> "l" <> ToString[#[[2]]] <> "m" <> ToString[#[[3]]] ->
          Join[ReIm[h @@ #], ReIm[Dr @@ #], 2]) & /@ modesLM]]];
Export[out, ds, "Datasets"];
Print["wrote ", out, " (", Length[ds] - 1, " datasets)"];
