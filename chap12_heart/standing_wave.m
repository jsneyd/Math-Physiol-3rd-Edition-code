%  -------------------------------------------------------------------
%
%   Symbolic calculations for the standing wave in a passive cardiac cable.
%
%   For Chapter 12, Section 12.4.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 

clear all
close all
clc

syms a1 a2 b L mu x l ratc rc re rgbyrc E;

a1 = mu - 1/E;
a2 = mu - E;
l = log(E)/L;

% Define the equations
ph = a1 * exp(l * x) + a2 * exp(-l * x);
phi = ratc * ph + b;
phe = -(1 - ratc) * ph + b;

eq1 = subs(phe,x,L) - mu * subs(phe,x,0);

% Solve for b
solb = solve(eq1, b);
pretty(simplify(solb))

phi = subs(phi,b,solb);
dphi = diff(phi,x);
dphiL = subs(dphi,x,L);
eq4 = subs(phi,x,L) - mu * subs(phi,x,0) + rgbyrc * dphiL;
pretty(simplify(eq4,150))   % This is the equation to solve for mu
