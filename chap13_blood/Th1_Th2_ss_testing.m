% this plots the Th1 Th2 steady state curves for Fig. 13.17
clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);


% set parameters from Table 13.5

sig = 5;
kap = 1;
rho = 1;
beta = 0.05;
alpha = 5;
g1 = 1;
g2 = 0.5;
mu = 3;
g = 2;
n = 4;

% select a value of S1 and S2
S1 = 0;
S2 = 0.2;

l = 3;
x1 = linspace(0,l,10000);
f1 = (-beta + mu*x1)./(alpha*x1.^n/(kap^n+x1.^n) +sig*S1/(rho+S1)) ;
x2=(g2./f1)-g2;

xx2 = linspace(0,l,10000);
f2 = (-beta + mu*xx2)./(alpha*xx2.^n/(kap^n+xx2.^n) +sig*S2/(rho+S2)) ;
xx1=(g1./f2)-g1;


figure(1)
plot(x1,x2,'--',xx1,xx2)
xlim([0,0.2])
ylim([0,l])


