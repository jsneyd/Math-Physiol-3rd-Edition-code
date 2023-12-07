%-------------------------------------------------------------------

% Matlab code to find analytic expressions for the flux through a
% saturating barrier model of an ion channel, with N binding sites, where
% every barrier is symmetric and identical.

% For Chapter 3, Exercise 3.6 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

clear all
close all
clc

syms x ci ce kon J k km konL konR 
syms v kbar

% number of binding sites. Note that with N binding sites there are 2
% endpoint equations, an additional N-1 equations for the binding sites
% and a conservation equation.

% 

N = 2;  
eq = sym('eq',[1 N-1]);
c = sym('c',[1 N]);

% Define the endpoint equations
eq0 = konL*ci*x - km*c(1) == J;
eqNN = k*c(N) - konR*ce*x == J;

% Define the internal equations
for i = 1:N-1
    eq(i) = k*c(i) - km*c(i+1) == J;
end

% Define the conservation equation
eqcons = x + sum(c) == 1;

% Solve for J in terms of ci and ce
sol = solve([eq0,eq,eqNN,eqcons],[x,c,J]);
pretty(simplify(sol.J))


% Make the rate constants functions of v = VF/RT
J  = subs(sol.J,[konL,konR,k,km],[kon*exp(v/(2*(N+1))),kon*exp(-v/(2*(N+1))),kbar*exp(v/(2*(N+1))),kbar*exp(-v/(2*(N+1)))]);
J = simplify(J)
pretty(J)

J_max = limit(J,ci,Inf)
pretty(J_max)

%% This bit is for calculating the derivative of Jmax

down = 1;
for i=2:N
    down = down + i*exp( (i-1)*v/(N+1) );
end

J_max = exp( (2*N-1)*v/(2*N+2) )/down;  % general expression for J_max

pretty(simplify(diff(J_max,v)))











