%-------------------------------------------------------------------

% Matlab code to find analytic expressions for the flux through a
% saturating barrier model of an ion channel, with 2 binding sites. Here,
% the channel can bind two ions at once.

% For Chapter 3, Section 3.4.3 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

clear all
close all
clc

syms k0 k1 k2 km1 km2 km3
syms p00 p01 p10 p11 ci ce

a1 = k0*ci + km3*ce;
a2 = km1 + k1 + km3*ce;
a3 = k2 + km2 + k0*ci;
a4 = k2+km1;

eq1 = -a1*p00 + km1*p10 + k2*p01 == 0;
eq2 = k0*ci*p00 - a2*p10 + km2*p01 + k2*p11 == 0;
eq3 = km3*ce*p00 + k1*p10 - a3*p01 + km1*p11 == 0;
eq4 = p00 + p01 + p10 + p11 - 1 == 0;

sols = solve([eq1,eq2 eq3,eq4],[p00,p01,p10,p11]);
p10 = sols.p10;
p01 = sols.p01;

J = p10*k1 - p01*km2;
[NJ,DJ] = numden(J);
pretty(simplify(NJ))
pretty(simplify(DJ))

simplify(limit(J,ci,Inf))
simplify(limit(J,ce,Inf))