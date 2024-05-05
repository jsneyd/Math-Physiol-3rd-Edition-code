%  -------------------------------------------------------------------
% 
%   Matlab code to make the plot Fig. 4.5 to determiine the eigenvalues as
%   the roots of Equation (4.66).
% 
%   For Chapter 4, Exercise 4.4 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  -------------------------------------------------------------------

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2);

L = 1;
gam = 1;
lam= [-1:.01:10];

f=tan(lam*L);
indx1 = find(lam<pi/2);
indx2 = find(lam>pi/2&lam<3*pi/2);
indx3 = find(lam>3*pi/2&lam<5*pi/2);
indx4 = find(lam>5*pi/2);

g = -lam/gam;

figure(1)
plot(lam(indx1),f(indx1),'r',lam(indx2),f(indx2),'r',lam(indx3),f(indx3),'r',lam(indx4),f(indx4),'r', ...
    lam,g,'b',lam,zeros(length(lam),1),'k--')
xlabel('\lambda')
axis([-1 10 -10 10])
box off
