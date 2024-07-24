%-------------------------------------------------------------------

% Matlab code to find analytic expressions for the flux through a
% saturating barrier model of an ion channel, with N binding sites

% For Chapter 3, Section 3.4.2 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

clear all
close all
clc

syms x ci ce k0 J
syms v kap0 

% N= number of binding sites. Note that with N binding sites there are 2
% endpoint equations, an additional N-1 equations for the binding sites
% and a conservation equation.

N = 3;  
eq = sym('eq',[1 N-1]);
k = sym('k',[1 N]);
km = sym('km',[1 N+1]);
kbar = sym('kbar',[1 N]);
kbarm = sym('kbarm',[1 N+1]);
c = sym('c',[1 N]);

% Define the endpoint equations
eq0 = k0*ci*x - km(1)*c(1) == J;
eqNN = k(N)*c(N) - km(N+1)*ce*x == J;

% Define the internal equations
for i = 1:N-1
    eq(i) = k(i)*c(i) - km(i+1)*c(i+1) == J;
end

% Define the conservation equation
eqcons = x + sum(c) == 1;

% Solve for J in terms of x, ci and ce
sol = solve([eq0,eq,eqNN],[c,J]);
pretty(simplify(sol.J))

% Include the conservation equation and solve for J in terms of ci and ce
sol = solve([eq0,eq,eqNN,eqcons],[x,c,J]);
pretty(simplify(sol.J))

% Make the rate constants functions of v = VF/RT
J  = subs(sol.J,[k0,k,km],[kap0*exp(v/(2*(N+1))),kbar*exp(v/(2*(N+1))),kbarm*exp(-v/(2*(N+1)))]);
pretty(simplify(J))

% For Exercise 3.5
pretty(simplify(limit(sol.J,ci,Inf)))

% For Exercise 3.6
pretty(simplify(limit(J,ci,Inf)))






