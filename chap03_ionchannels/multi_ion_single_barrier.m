%-------------------------------------------------------------------

% Matlab code to find analytic expressions for the flux through a
% saturating barrier model of an ion channel, with N binding sites. Here,
% the channel can bind two ions at once.

% For Chapter 3, Section 1.4.3 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

clear all
close all
clc

syms k12 k13 k14 k21 k23 k24 k31 k32 k34 k41 k42 k43 kf kon
syms p1 p2 p3 p4 ci ce

a1 = 2*kf;
a2 = kon*ce + kf + k42;
a3 = kon*ci + kon*ce;
a4 = kon*ci + k24 + kf;

eq1 = -a1*p1 + kon*ce*p2 + kon*ci*p4 == 0;
eq2 = kf*p1 - a2*p2 + kon*ci*p3 + k24*p4 == 0;
eq3 = kf*p2 - a3*p3 + kf*p4 == 0;
eq4 = p1 + p2 + p3 + p4 - 1 == 0;

sols = solve([eq1,eq2 eq3,eq4],[p1,p2,p3,p4]);
p2 = sols.p2;
p4 = sols.p4;

J = p2*k42 - p4*k24;
[NJ,DJ] = numden(J);
pretty(simplify(NJ))

limit(J,kf,Inf)
limit(J,ci,Inf)
limit(J,ce,Inf)