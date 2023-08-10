clear all
close all
clc

syms B A N M r k sigma rN aN dN aM rA dA m rM aM dM K

A = r/k;
F = B;

sigma = 0;
rN = 0;
rM = 0;
%rA = 0;

eq1 = sigma + rN*N == aN*F*N + dN*N;
eq2 = aN*F*N + aM*F*M + rA*F*A == dA*A + m*(1-F)*A;
eq3 = m*(1-F)*A + rM*M == aM*F*M + dM*M;

sols = solve([eq1 eq2 eq3],[B N M]);

pretty(simplify(sols.B,150))
pretty(simplify(sols.M,150))

