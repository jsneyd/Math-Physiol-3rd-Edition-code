

clear all
close all
clc

syms G V0 V1 V2

V = V0 + G*V1;

wun = 1 - 2*G*V*(1+V/G);
cc = collect(simplify(wun),G);
cf = coeffs(cc,G);
cf1 = simplify(cf(1))
V0sol = solve(cf1==0,V0)
cf2 = simplify(subs(cf(2),V0,V0sol))
V1sol = solve(cf2==0,V1)

