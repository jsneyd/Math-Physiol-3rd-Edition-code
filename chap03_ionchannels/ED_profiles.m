%  -------------------------------------------------------------------
%
%   Plot profiles from the PNP equations.
%
%   For Chapter 3, Section 3.3.1 of
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


y = linspace(0,1,100);
ui = 50/550;
ue = 500/550;
v = 1;
J1 = v*(ui - ue*exp(-v))/(1-exp(-v));
K1 = ui - J1/v;

cprofile = J1/v + K1*exp(v*y);
vprofile= v-v*y;
profilex = y;
cprofilelong = ui +(ue-ui)*y;
vprofilelong = -(v/log(ue/ui))*log(ui/ue + (1-ui/ue)*y);

figure(1)
    plot(y,cprofile,y,vprofile)
    xlabel('y')
    legend('u_1','\psi')

figure(2)
    plot(y,cprofilelong,y,vprofilelong)
    xlabel('y')
    legend('u_1','\psi')