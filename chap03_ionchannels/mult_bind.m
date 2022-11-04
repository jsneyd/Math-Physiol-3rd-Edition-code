clear all
close all
clc

syms k0 k1 k2 k3 km1 km2 km3 km4
syms c0 c1 c2 c3 c4 x

% n=3 case
eq1 = k0*c0*x - km1*c1 ==  k1*c1 - km2*c2;
eq2 = k1*c1 - km2*c2 == k2*c2 - km3*c3*x;
cons = x+c1+c2==1;

[c1,c2,x] = solve([eq1,eq2,cons],[c1,c2,x]);
J = simplify(k0*c0*x - km1*c1,'Steps',50);
pretty(J)

% n=4 case
clear all
syms k0 k1 k2 k3 km1 km2 km3 km4 K1 K2 K3
syms c0 c1 c2 c3 c4 x

eq1 = k0*c0*x - km1*c1 ==  k1*c1 - km2*c2;
eq2 = k1*c1 - km2*c2 == k2*c2 - km3*c3;
eq3 = k2*c2 - km3*c3 == k3*c3 - km4*c4*x;
cons = x+c1+c2+c3==1;

[c1,c2,c3,x] = solve([eq1,eq2,eq3,cons],[c1,c2,c3,x]);
J = simplify(k0*c0*x - km1*c1,'Steps',50);
pretty(J)