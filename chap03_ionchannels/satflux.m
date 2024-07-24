%  -------------------------------------------------------------------
%
%   Flux in a saturating one-ion model.
%
%   For Chapter 3, Section 3.4.2 of
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
   'defaultlinelinewidth', 2.0);

k0 = 1;
k1 = 1;
km1 = 1;
km2 = 1;

alpha1 = (km1+k1)/(k0*k1);
beta1 = 1/k1;
gamma1 = km2/(k0*k1);

ce = 1;

ci = linspace(0,10,100);
J = ( ci - ce*(km1*km2/(k0*k1)) )./(alpha1 + beta1*ci + gamma1*ce);

plot(ci,J)
xlabel('c_i')
ylabel('J')
