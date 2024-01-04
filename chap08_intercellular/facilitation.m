%  -------------------------------------------------------------------
%
%   Program to compute facilitation in a model of the residual bound
%   calcium hypothesis of synaptic plasticity.
%
%   For Chapter 8, Fig. 8.9 of
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
tpcp=200;
k1=3.75e-3;
k2=2.5e-3;
k3=5e-4;

km1=4e-4;  %time units are ms
km2=1e-3;
km3=0.1;
K1=km1/k1;
K2 = km2/k2;
K3 = km3/k3;

T=[1/0.186:0.001:10000];
a1=exp(-km1*(T+tpcp/K1));  
a2=exp(-km2*(T+tpcp/K2));
a3=exp(-km3*(T+tpcp/K3));

Fmax = (1./(1-a1)).*(1./(1-a2)).*(1./(1-a3)); 
    
semilogx(1000./T,Fmax)% the factor 1000 is to convert to Hertz  = s^{-1}
xlabel('Stimulus Frequency (Hz)')
ylabel('Facilitation')
axis([0.1 100 0 7])