(* ::Package:: *)

(* Expansion about the pole (alpha = 0) of the directly regularised SS source at r = r0, for all r0:
       d2R[hS,hS] - Box[hSS] - GD[hSS] = c0(beta, r0) + c2(beta, r0) alpha^2 + c4(beta, r0) alpha^4 + ...
   (Box: plain Lorenz operator; GD: its gauge-damping term; Box + GD = the gauge-damped operator E0).

   At r = r0 the components of d2R and Box[hSS] diverge like rho^-4 at the pole while the difference
   is finite; evaluating the closed forms there needs ever more precision. With Delta r = 0 exactly,
   alpha -> 0 is an ordinary limit, so the difference is expanded in alpha; d2R-SS-quadrature.wl uses
   the expansion for the r0 row below a small alpha.

   Method (block by block, nothing inlined): the module bodies of d2R[mu,nu] (d2R.wl), BoxhSS[mu,nu] and
   GDhSS[mu,nu] (hSS.wl) are re-evaluated with r = r0 = R0 symbolic (Delta r = 0), the surds as symbols
   (rho = 2 q y Sin[alpha/2], q = r0 Sqrt[(r0-2)/(r0-3)]; chi = y^2; Omega0 = w = r0^(-3/2)), and each
   small block (the regularised hbar / hSS and their derivatives, 4x4 each) converted to a truncated
   series in alpha as soon as it is built. The modules' own closed-form combinations then combine
   series (truncated at every product). The coefficients are reduced (q^2, y^2, w^2 -> rational) and
   simplified at the end; the negative powers of alpha must cancel identically.

   Usage: wolframscript -file d2R-minus-box-hSS-r0-series.wl [components] [blockOrder] [nkernels]
     components: e.g. "{1}" (default all 10); blockOrder: order in alpha of the blocks (default 14)
   Output: d2R-minus-box-hSS-r0-series.m: {"components" -> ..., "order" -> n,
     "coefficients" -> {{c0, c2, ...} per component}} in R0 and b (beta). *)

$here = DirectoryName[$InputFileName];
args = Rest[$ScriptCommandLine];
comps = If[Length[args] >= 1, ToExpression[args[[1]]], Range[10]];
blockOrder = If[Length[args] >= 2, ToExpression[args[[2]]], 14];
nk = If[Length[args] >= 3, ToExpression[args[[3]]], Min[10, Length[comps]]];
compNames = {"tt", "tr", "rr", "ta", "tb", "ra", "rb", "aa", "ab", "bb"};

(* the source text of a module "name[r_,r0_]:=" ... "]" *)
moduleLines[file_, name_] := Module[{lines = StringSplit[ReadString[FileNameJoin[{$here, file}]], "\n"], k0, k1},
   k0 = First@FirstPosition[lines, s_String /; StringStartsQ[s, name <> "[r_,r0_]:="]];
   k1 = k0 + First@FirstPosition[lines[[k0 + 1 ;;]], s_String /; s === "]"];
   lines[[k0 ;; k1]]];

(* series version: Delta r = 0, surds as symbols, blocks -> series right after they are built, the series
   branches / switches / Compile replaced by the closed-form combinations *)
seriesModule[file_, name_, blockLHS_, numLHS_, result_] := Module[{l = moduleLines[file, name], out = {}, skip = False},
   Do[Which[
      StringStartsQ[s, name <> "[r_,r0_]:="], AppendTo[out, name <> "Ser[r_,r0_,nOrd_]:="],
      StringStartsQ[s, "\\[Chi]="], AppendTo[out, "\\[Chi]=y^2;"],
      StringStartsQ[s, "\\[Rho]="], AppendTo[out, "\\[Rho]=2 q y Sin[\\[Alpha]/2];"],
      StringStartsQ[s, "\\[CapitalOmega]0="], AppendTo[out, "\\[CapitalOmega]0=w;"],
      StringStartsQ[s, "{" <> blockLHS <> ","], AppendTo[out, s];
        With[{names = "{" <> First[StringSplit[StringDrop[s, 1], "}="]] <> "}"},
          AppendTo[out, names <> "=Map[Series[#,{\\[Alpha],0,nOrd}]&," <> names <> ",{3}];"]],
      (* alpha and beta global (not Module-local), so that the series and coefficients use them *)
      StringContainsQ[s, "\\[Alpha],\\[Beta],"] && !StringContainsQ[s, "="], AppendTo[out, StringReplace[s, "\\[Alpha],\\[Beta]," -> ""]],
      StringContainsQ[s, "series}="] || StringContainsQ[s, "series,"] && StringContainsQ[s, "}={"], Null,   (* series branches *)
      StringContainsQ[s, "=If[\\[Alpha]<\\[Epsilon]\\[Alpha],"] || StringStartsQ[s, "\\[Epsilon]\\[Alpha]="], Null,
      StringStartsQ[s, "{{r,r0},{Compile["], skip = True; AppendTo[out, result],
      skip && s =!= "]", Null,
      True, AppendTo[out, s]], {s, l}];
   ToExpression[StringRiffle[out, "\n"]]];

d2RStart = AbsoluteTime[];
seriesModule["d2R.wl", "d2R\\[Mu]\\[Nu]", "hbarreg", "",
  "{d2R\\[ScriptCapitalA]ttnum+d2R\\[ScriptCapitalB]ttnum,d2R\\[ScriptCapitalA]trnum+d2R\\[ScriptCapitalB]trnum,d2R\\[ScriptCapitalA]rrnum+d2R\\[ScriptCapitalB]rrnum,d2R\\[ScriptCapitalA]t\\[Theta]num+d2R\\[ScriptCapitalB]t\\[Theta]num,d2R\\[ScriptCapitalA]t\\[Phi]num+d2R\\[ScriptCapitalB]t\\[Phi]num,d2R\\[ScriptCapitalA]r\\[Theta]num+d2R\\[ScriptCapitalB]r\\[Theta]num,d2R\\[ScriptCapitalA]r\\[Phi]num+d2R\\[ScriptCapitalB]r\\[Phi]num,d2R\\[ScriptCapitalA]\\[Theta]\\[Theta]num+d2R\\[ScriptCapitalB]\\[Theta]\\[Theta]num,d2R\\[ScriptCapitalA]\\[Theta]\\[Phi]num+d2R\\[ScriptCapitalB]\\[Theta]\\[Phi]num,d2R\\[ScriptCapitalA]\\[Phi]\\[Phi]num+d2R\\[ScriptCapitalB]\\[Phi]\\[Phi]num}"];
boxResult = "{Boxhttnum,Boxhtrnum,Boxhrrnum,Boxht\\[Theta]num,Boxht\\[Phi]num,Boxhr\\[Theta]num,Boxhr\\[Phi]num,Boxh\\[Theta]\\[Theta]num,Boxh\\[Theta]\\[Phi]num,Boxh\\[Phi]\\[Phi]num}";
seriesModule["hSS.wl", "BoxhSS\\[Mu]\\[Nu]", "hreg", "", boxResult];
seriesModule["hSS.wl", "GDhSS\\[Mu]\\[Nu]", "hreg", "", boxResult];
Print["series modules defined"];

(* The blocks contain Log[rho] (rho = 2 q y Sin[alpha/2]), so the series coefficients carry Log[c alpha^k]
   factors whose arguments differ in form between d2R and Box/GD; split each into k LA + Log[c]
   (LA = Log[alpha]). The coefficients of the negative powers of alpha must vanish; the difference is
   then c(0,0) + c(0,1) Log[alpha] + O(alpha): finite apart from a weak (integrable) log at the pole. *)
splitLogs[e_] := e /. Log[x_] /; !FreeQ[x, \[Alpha]] :> Module[{k = Exponent[x, \[Alpha]], rest},
     rest = Simplify[x/\[Alpha]^k];
     If[FreeQ[rest, \[Alpha]], k LA + Log[rest], Print["  WARNING: log argument not of the form c alpha^k: ", Short[x]]; Log[x]]];
coefficientTable[s_, pmin_, pmax_] := Module[{e = Expand[splitLogs[Normal[s]]]},
   Table[{p, j, Coefficient[Coefficient[e, \[Alpha], p], LA, j]}, {p, pmin, pmax}, {j, 0, 2}]];

(* Diagnostic: SERIES_DIAG = r0 (e.g. 179/10): everything but alpha numeric (50 digits; beta = 6/5), all
   components; checks the cancellation and compares with the closed forms at 100 digits. Fast. *)
If[StringQ[Environment["SERIES_DIAG"]],
  Module[{R0n = ToExpression[Environment["SERIES_DIAG"]], bn = 6/5, prec = 50, t0 = AbsoluteTime[], sers, diffs, ok = True,
     fd, fb, fg, body},
    q = SetPrecision[R0n Sqrt[(R0n - 2)/(R0n - 3)], prec]; y = SetPrecision[Sqrt[1 - Sin[bn]^2/(R0n - 2)], prec];
    w = SetPrecision[R0n^(-3/2), prec]; \[Beta] = SetPrecision[bn, prec];
    sers = With[{o = blockOrder}, {d2R\[Mu]\[Nu]Ser[R0n, R0n, o], BoxhSS\[Mu]\[Nu]Ser[R0n, R0n, o], GDhSS\[Mu]\[Nu]Ser[R0n, R0n, o]}];
    Print["numeric modules evaluated as series in ", Round[AbsoluteTime[] - t0], " s"];
    diffs = sers[[1]] - sers[[2]] - sers[[3]];
    Get[FileNameJoin[{$here, "d2R.wl"}]]; Get[FileNameJoin[{$here, "hSS.wl"}]];
    body[cf_] := With[{fn = FirstCase[List @@ cf, _Function]}, fn[#1, #2] /. If[_, _, e_] :> e &];
    fd = body /@ d2R\[Mu]\[Nu][R0n, R0n][[2]]; fb = body /@ BoxhSS\[Mu]\[Nu][R0n, R0n][[2]]; fg = body /@ GDhSS\[Mu]\[Nu][R0n, R0n][[2]];
    Do[Module[{s = diffs[[k]], ord, tab, d2tab, bad, scale, sr, cf, rel},
       ord = If[Head[s] === SeriesData, s[[5]]/s[[6]], Infinity];
       tab = Flatten[coefficientTable[s, -6, ord - 1], 1];
       d2tab = Flatten[coefficientTable[sers[[1, k]], -6, -1], 1];
       scale = Max[Abs[N[d2tab[[All, 3]]]], 10^-300];
       bad = Select[tab, #[[1]] < 0 &];
       Print[compNames[[k]], ": order alpha^", ord, "; largest d2R coefficient of a negative power ", ScientificForm[scale, 3],
         "; largest leftover negative power ", ScientificForm[Max[Abs[N[bad[[All, 3]]]], 0], 3],
         "; c(0,0) = ", N[Select[tab, #[[;; 2]] == {0, 0} &][[1, 3]], 12], ", c(0,1) [Log alpha] = ", N[Select[tab, #[[;; 2]] == {0, 1} &][[1, 3]], 12]];
       If[Max[Abs[N[bad[[All, 3]]]], 0] > 10^-30 scale, ok = False];
       Do[
         sr = Total[(#[[3]] al^#[[1]] Log[al]^#[[2]]) & /@ Select[tab, #[[1]] >= 0 &]];
         cf = N[fd[[k]][SetPrecision[al, 100], SetPrecision[bn, 100]] - fb[[k]][SetPrecision[al, 100], SetPrecision[bn, 100]] -
              fg[[k]][SetPrecision[al, 100], SetPrecision[bn, 100]], 30];
         rel = Abs[sr/cf - 1];
         Print["   alpha ", N[al], ": closed form ", N[cf, 16], "  series ", N[sr, 16], "  rel. diff ", ScientificForm[N[rel], 2]];
         If[!TrueQ[al == 1/10000 && rel < 10^-12] && al == 1/10000, ok = False], {al, {1/10000, 1/1000, 1/100}}]], {k, comps}];
    Print[If[ok, "DIAGNOSTIC PASSED", "DIAGNOSTIC FAILED"], " in ", Round[AbsoluteTime[] - t0], " s"];
    Quit[If[ok, 0, 1]]]];

rel = {q2 -> R0^2 (R0 - 2)/(R0 - 3), y2 -> 1 - Sin[b]^2/(R0 - 2), w2 -> R0^-3};
reduce[e_] := Expand[e] /. {q^n_Integer :> q^Mod[n, 2] q2^Quotient[n, 2], y^n_Integer :> y^Mod[n, 2] y2^Quotient[n, 2],
     w^n_Integer :> w^Mod[n, 2] w2^Quotient[n, 2]} /. rel;
restore = {q -> R0 Sqrt[(R0 - 2)/(R0 - 3)], y -> Sqrt[1 - Sin[b]^2/(R0 - 2)], w -> R0^(-3/2)};

LaunchKernels[nk];
ParallelEvaluate[$here = #] &[$here];
DistributeDefinitions["Global`"];
t0 = AbsoluteTime[];
(* the modules evaluate all 10 components at once; the reduction is per component *)
ser = With[{o = blockOrder}, {d2R\[Mu]\[Nu]Ser[R0, R0, o], BoxhSS\[Mu]\[Nu]Ser[R0, R0, o], GDhSS\[Mu]\[Nu]Ser[R0, R0, o]}] /. \[Beta] -> b;
Print["modules evaluated as series in ", Round[AbsoluteTime[] - t0], " s"];
diff = ser[[1]] - ser[[2]] - ser[[3]];
result = ParallelTable[
   Module[{s = diff[[k]], pw, coef, ord},
     ord = If[Head[s] === SeriesData, s[[5]]/s[[6]], Infinity];
     s = Expand[splitLogs[Normal[s]]];
     pw = If[s === 0, {0}, Union[Exponent[#, \[Alpha]] & /@ If[Head[s] === Plus, List @@ s, {s}]]];
     coef = Association[Table[{p, j} -> Together[reduce[Coefficient[Coefficient[s, \[Alpha], p], LA, j]]], {p, pw}, {j, 0, 2}]];
     {k, ord, pw, AllTrue[Keys[coef], #[[1]] >= 0 || coef[#] === 0 &],
      AllTrue[Select[Keys[coef], OddQ[#[[1]]] &], coef[#] === 0 &],
      (* {p, j} -> coefficient of alpha^p Log[alpha]^j, 0 <= p < order *)
      Association[KeyValueMap[#1 -> Simplify[#2 /. restore] &, KeySelect[DeleteCases[coef, 0], #[[1]] >= 0 && #[[1]] < ord &]]]}],
   {k, comps}, Method -> "FinestGrained"];
Do[Print["component ", compNames[[r[[1]]]], ": truncation order alpha^", r[[2]], ", powers ", r[[3]],
    "; negative powers vanish: ", r[[4]], "; odd powers vanish: ", r[[5]], "; coefficients kept: ", Length[r[[6]]]], {r, result}];
Print["reduced in ", Round[AbsoluteTime[] - t0], " s"];

(* validation at r0 = 17.9, beta = 6/5 against the closed forms at 100 digits (d2R.wl d2R, hSS.wl Box, GD) *)
Get[FileNameJoin[{$here, "d2R.wl"}]]; Get[FileNameJoin[{$here, "hSS.wl"}]];
body[cf_] := With[{fn = FirstCase[List @@ cf, _Function]}, fn[#1, #2] /. If[_, _, e_] :> e &];
fd = body /@ d2R\[Mu]\[Nu][179/10, 179/10][[2]]; fb = body /@ BoxhSS\[Mu]\[Nu][179/10, 179/10][[2]]; fg = body /@ GDhSS\[Mu]\[Nu][179/10, 179/10][[2]];
ok = True;
Do[With[{k = r[[1]], cs = r[[6]]},
   Do[Module[{x = SetPrecision[al, 100], yb = SetPrecision[6/5, 100], cf, sr},
      cf = N[fd[[k]][x, yb] - fb[[k]][x, yb] - fg[[k]][x, yb], 30];
      sr = N[Total[KeyValueMap[#2 a^#1[[1]] Log[a]^#1[[2]] &, cs]] /. {R0 -> 179/10, b -> 6/5, a -> al}, 30];
      Print["  ", compNames[[k]], " alpha ", N[al], ": closed form ", N[cf, 16], "  series ", N[sr, 16], "  rel. diff ", ScientificForm[N[Abs[sr/cf - 1]], 2]];
      If[al == 1/10000 && !TrueQ[Abs[sr/cf - 1] < 10^-12], ok = False]], {al, {1/10000, 1/1000, 1/100}}]], {r, result}];
If[ok && AllTrue[result, #[[4]] &],
  Put[{"components" -> compNames[[comps]], "order" -> (result[[All, 2]]), "coefficients" -> result[[All, 6]], "form" -> "{p, j} -> c: c alpha^p Log[alpha]^j"},
    FileNameJoin[{$here, "d2R-minus-box-hSS-r0-series.m"}]];
  Print["wrote d2R-minus-box-hSS-r0-series.m in ", Round[AbsoluteTime[] - t0], " s"],
  Print["NOT written (validation ok: ", ok, ")"]];
