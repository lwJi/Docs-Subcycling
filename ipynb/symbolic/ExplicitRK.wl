(* ::Package:: *)

(* ::Section:: *)

(*ExplicitRK*)

(* ::Input::Initialization:: *)

Y[s_] :=
    y0 + h Sum[Subscript[a, s, i] f[Y[i], t0 + h Subscript[c, i]], {i,
         1, s - 1}];

y1[s_] :=
    y0 + h Sum[Subscript[b, i] f[Y[i], t0 + h Subscript[c, i]], {i, 1,
         s}];

k[i_] :=
    f[Y[i], t0 + h Subscript[c, i]];

Protect[y0, t0, h, a, b, c];

Subscript[c, 1] = 0;

(* ::Section:: *)

(*Rules*)

(* ::Input::Initialization:: *)

autonomousRules1 = {Subscript[a, 2, 1] -> Subscript[c, 2], Subscript[
    a, 3, 1] + Subscript[a, 3, 2] -> Subscript[c, 3], Subscript[a, 4, 1] 
    + Subscript[a, 4, 2] + Subscript[a, 4, 3] -> Subscript[c, 4], Derivative[
    i_, 1][f][y0, t0] -> 0, Derivative[i_, 2][f][y0, t0] -> 0, Derivative[
    i_, 3][f][y0, t0] -> 0};

(* ::Input::Initialization:: *)

regularRK4Rules = {Subscript[a, 2, 1] -> 1/2, Subscript[a, 3, 1] -> 0,
     Subscript[a, 3, 2] -> 1/2, Subscript[a, 4, 1] -> 0, Subscript[a, 4, 
    2] -> 0, Subscript[a, 4, 3] -> 1, Subscript[c, 2] -> 1/2, Subscript[c,
     3] -> 1/2, Subscript[c, 4] -> 1, Subscript[b, 1] -> 1/6, Subscript[b,
     2] -> 1/3, Subscript[b, 3] -> 1/3, Subscript[b, 4] -> 1/6};

(* ::Input::Initialization:: *)

autonomousRules2 =
    Module[{dydt1Rules, dydt2Rules, dydt3Rules},
        dydt1Rules = {f[y0_, t0_] -> y'[t0]};
        dydt2Rules = {Derivative[1, 0][f][y0_, t0_] -> y''[t0] / f[y0,
             t0]};
        dydt3Rules = {Derivative[2, 0][f][y0_, t0_] -> (y'''[t0] - f[
            y0, t0] Derivative[1, 0][f][y0, t0] ^ 2) / f[y0, t0] ^ 2};
        Join[dydt1Rules, dydt2Rules, dydt3Rules /. dydt2Rules // Simplify
            ]
    ];
