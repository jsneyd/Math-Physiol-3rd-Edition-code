
%-------------------------------------------------------------------
%
% Discrete release standing waves.
%
% For Chapter 7, Section 7.8 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
% Written by James Keener and James Sneyd
%
%-------------------------------------------------------------------

clear all
close all
clc

set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0);



Abyk =5;
beta = [0.01:.1:5];

astar = 1-sqrt(1/(1+Abyk));
Af = Abyk.*beta.*sinh(beta);
l0 = cosh(beta);
lf = l0+Af/2;
mu0 = l0-sqrt(l0.^2-1);
muf = lf-sqrt(lf.^2-1);

C = Af./(Af+2*l0-2);

a = C.*(muf-1)./(muf-1./mu0);
b = C.*(1-1./mu0)./(muf-1./mu0);
bf = C-b.*muf;

figure(1)
    plot(beta,a,'r',beta,bf,'b','linewidth',2)
    xlabel('\beta','fontsize',18)
    ylabel('c^*/c_e','fontsize',18)
    text(2,.4,'standing waves','fontsize',20)
    text(.2, .1,'traveling waves','fontsize',20)

%output=[beta' a' bf'];
%save discrete.dat output -ascii
