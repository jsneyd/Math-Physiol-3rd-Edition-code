function assay

clear all
close all
clc

syms R V0 V1 V2

V = V0 + R*V1;
ff = taylor(1/V, R, 'Order', 4);
ee = taylor(1-exp(-R/(2*V)),R,'Order',4);

wun = 1 - 2*V/R*ee*(1+R*V);
cc = collect(simplify(wun),R);
cf = coeffs(cc,R);
cf1 = simplify(cf(1))
V0sol = solve(cf1==0,V0)
cf2 = simplify(subs(cf(2),V0,V0sol))
V1sol = solve(cf2==0,V1)

end