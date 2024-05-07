%    -------------------------------------------------------------------
%
%     Plot the function that controls light adaptation in the model of LA
%     in cones. Compare to data.
%
%     For Chapter 19, Section 19.2.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);

vstar = 35.7;
s1=1.59/vstar;
s2=1130;
vK=-13/vstar ; 
tauy=0.07; 
k1=35.4 ;
gam=303 ;
delta=5 ;
kappa=0.1;
eta=52.5 ;
tau1=0.012 ;
taum=0.016 ;
tauz=0.04;

y = linspace(0.2,1,100);

phi=(y.*exp(vK*(1-y))).^(1/3).*(delta+(gam-delta)*eta*(exp(-vK*(1-y)/s1)-1)./ ...
        (s2*k1+eta*(exp(-vK*(1-y)/s1)-1)));
A=4+84./(1+(y/0.34).^4); 

plot(y,phi,'r',y,A,'--b')
xlabel('y')
ylabel('[cGMP]_{dark} s^{-1}')
legend('Model','Data')