%-------------------------------------------------------------------

% Matlab code to find analytic expressions for the flux through a
% nonsaturating barrier model of an ion channel, with N binding sites

% For Chapter 3, Section 1.4.1 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

clear all
close all
clc

syms x ci ce k0 M 

% N is the number of binding sites. Note that with N binding sites you have 2
% endpoint equations and an additional N-1 equations for the binding sites.

% The text is written in terms of N barriers, and thus N-1 binding sites,
% but it's slightly easier to write the code for N binding sites. 

N = 3;  
eq = sym('eq',[1 N-1]);
k = sym('k',[1 N]);
km = sym('km',[1 N+1]);
c = sym('c',[1 N]);

% Define the endpoint equations
eq0 = k0*ci - km(1)*c(1) == M;
eqNN = k(N)*c(N) - km(N+1)*ce == M;

% Define the internal equations
for i = 1:N-1
    eq(i) = k(i)*c(i) - km(i+1)*c(i+1) == M;
end

% Solve for J
sol = solve([eq0,eq,eqNN],[c,M]);
pretty(simplify(sol.M))
