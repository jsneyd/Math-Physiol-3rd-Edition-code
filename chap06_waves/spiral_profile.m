% shooting to find spiral wave profile

 
%
%   For Chapter 6, Figure 1.14,   of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 

 clear all; close all; clc;
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);  
 
global c w 



[r,ph] = ode23(@oderhs,...)
    ps = ph./(sqrt(r.^2-ph.^2));  
     
function out = derhs(r,ph)
global c w 
dph = r*(c-w*sqrt(r^2-ph^2));
end
