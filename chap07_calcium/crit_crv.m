
%  -------------------------------------------------------------------
%
%   Plot of formula 8.153 for figure 8.28
%
%   For Chapter 8, 
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

Kb = [0.01:.01:10];
ac = -((3*Kb.^2+2*Kb).*log((Kb+1)./Kb)-(3*Kb+1/2))./((Kb+1/2).*log((Kb+1)./Kb)-1)/2;

semilogx(Kb,ac)
xlabel('K_b')
ylabel('a_{c,min}')
