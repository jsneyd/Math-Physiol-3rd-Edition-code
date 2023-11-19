%-------------------------------------------------------------------

% Matlab code for determining the flux through a Na/Ca exchanger.

% For Chapter 2, Section 2.3.3 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

close all
clear all
clc

syms x1 x2 y1 y2 k1 k2 k3 k4 km1 km2 km3 km4
syms ni ne ci ce

% Define the system of differential equations
eq1 = km1*ni^3*x2 + k4*y1 - (k1*ci + km4)*x1 == 0;
eq2 = km2*y2 + ci*k1*x1 - (km1*ni^3 + k2)*x2 == 0;
eq3 = km4*x1 + k3*ne^3*y2 - (k4 + km3*ce)*y1 == 0;
eq4 = x1 + x2 + y1 + y2 - 1 == 0;

% Solve the system of equations
sol = solve([eq1, eq2, eq3, eq4], [x1, x2, y1, y2]);

% Assign the solutions to variables
x1 = sol.x1;
x2 = sol.x2;
y1 = sol.y1;
y2 = sol.y2;

% Define the flux J
J = k4*y1 - km4*x1;

[N,D] = numden(J);
pretty(simplify(N,150))
pretty(simplify(D,150))
