(* ::Package:: *)

(* Two-stage, second-order RK *)


Y[s_]:=Subscript[y, n]+h Sum[a[s,i]f[Y[i]],{i,1,s-1}];


Y[1]


Y[2]


Series[f[Y[2]], {Subscript[y, n], 0, 1}]


dtRules = {
Derivative[1][y][x_]->f[y],
Derivative[2][y][x_]->f[y]f'[y],
Derivative[3][y][x_]->f[y]f'[y]f'[y]+f[y]f[y]f''[y]
}


Series[y[Subscript[t, n]+h], {h, 0, 2}]


Series[y[Subscript[t, n]+h], {h, 0, 2}]/.dtRules


D[f[x],x]==f'[x]
