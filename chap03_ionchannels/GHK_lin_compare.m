
%  -------------------------------------------------------------------
%
%   Comparing the GHK and linear IV curves.
%
%   For Chapter 3, Section 3.1.1 of
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

FRT = 1/25.8;
nai = 50;
nae = 437;
vna = (1/FRT)*log(nae/nai);
gna = 0.01;
pnaF = -gna*vna/(nai-nae);      %  to guarantee the curves meet at V=0

ki =  397;
ke = 20;
gk   = 0.367;
vk = (1/FRT)*log(ke/ki);
pkF   = -gk*vk/(ki-ke);         % this is fixed at value determined by normal kae

x = linspace(-40,100,100);                        
linIV = gna*(x-vna);
GHKIV = pnaF*FRT*x.*(nai - nae*exp(-x*FRT))./(1-exp(-x*FRT));

x1 = linspace(10,50,100); 
linVr = (gna* vna + gk* (1/FRT)*log(x1/ki))/(gna + gk);
GHKVr = (1/FRT) * log( (pnaF*nae + pkF*x1)/(pnaF*nai + pkF*ki) );   % Note that F cancels out

figure(1)
    plot(x,linIV,x,GHKIV)
    legend('linear','GHK')
    xlabel('V (mV)')
    ylabel('I_{Na}')

figure(2)
    plot(x1,linVr,x1,GHKVr)
    legend('linear','GHK')
    xlabel('[K^+]_e (mM)')
    ylabel('V_r (mV)')