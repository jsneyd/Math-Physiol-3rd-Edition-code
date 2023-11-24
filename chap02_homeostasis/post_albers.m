%-------------------------------------------------------------------

% Matlab code to compute the flux through the Post-Albers 
% model of the Na/K ATPase.

% For Chapter 2, Section 2.4.3 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

clear all
close all
clc

syms x1 x2 x3 y1 y2 y3 z1 z2
syms k1 k2 k3 k4 k5 k6 k7 k8
syms K1 K2 K3 K4 K5 K6 K7 K8
syms ADP ATP Pi Ni Ne Ki Ke

% Often it's easier to work with the k1 and K1, instead of the km1, etc.
km1 = K1*k1;
km2 = K2*k2;
km3 = K3*k3;
km4 = K4*k4;
km5 = K5*k5;
km6 = K6*k6;
km7 = K7*k7;
km8 = K8*k8;

eq1 = km1*Ki^2*x2 + k8*ATP*z1 - (km8+k1)*x1;        % x1 equation
eq2 = km2*x3 + k1*x1 - (km1*Ki^2 + k2*Ni^3)*x2;     % x2 equation
eq3 = km3*z2*ADP + k2*Ni^3*x2 - (k3+km2)*x3;        % x3 equation
eq4 = km4*y3 + k3*x3 - (k4 + km3*ADP)*z2;           % z2 equation
eq5 = k4*z2 + km5*Ne^3*y2 - (k5 + km4)*y3;          % y3 equation
eq6 = k5*y3 +km6*y1 - (km5*Ne^3+k6*Ke^2)*y2;        % y2 equation
eq7 = k6*Ke^2*y2 + km7*Pi*z1 - (km6 + k7)*y1;       % y1 equation
eq8 = x1 + x2 + x3 + y1 + y2 + y3 + z1 + z2 - 1;

sols = solve([eq1==0,eq2==0,eq3==0,eq4==0,eq5==0,eq6==0,eq7==0,eq8==0],[x1,x2,x3,y1,y2,y3,z1,z2]);
z2 = sols.z2;
y3 = sols.y3;

J = k4*z2 - km4*y3;
[N,D] = numden(J);
pretty(simplify(N))

