close all
clear all
clc
 
syms k kp km K Kd 
syms ce pe ci pi se si C0

kp = k/Kd;
km = kp*K;

eq1 = k*pe - k*pi + kp*si*ci - km*pi;
eq2 = k*pi - k*pe + kp*se*ce - km*pe;
eq3 = k*ce - k*ci + km*pi - kp*si*ci;
eq4 = k*ci - k*ce + km*pe - kp*se*ce;
eq5 = pi+pe+ci+ce-C0;

[ci,ce,pi,pe] = solve([eq2,eq3,eq4,eq5],[ci,ce,pi,pe])

J = simplify(km*pi - kp*si*ci,'Steps',150);
[up,down] = numden(J);
pretty(simplify(up))
pretty(simplify(down))

simplify(diff(J,si));
