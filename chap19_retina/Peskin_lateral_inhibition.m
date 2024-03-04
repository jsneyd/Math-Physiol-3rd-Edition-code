
% The Matlab/Octave file used to generate the image in Fig. 15 of Chapter 19
% Solutions of the Peskin lateral inhibition  model 
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);

%parameters
lam = 1;
k = 20;

t = [0:.01:3];

R = 1/(lam+1) -(k-1)/(k-lam-1)*exp(-k*t) + lam*k/((k-lam-1)*(lam+1))*exp(-(lam+1)*t);  

x = [-2:.01:2];
E = (x>=0);
sql = sqrt(lam+1);
Ix=lam/(lam+1)*((x>=0).*(1-exp(-x*sql)/2) +(x<0).*exp(x*sql)/2);


figure(1)
plot(t,R)
box off
xlabel('t')
ylabel('R(t)')

figure(2)
plot(x,E,'--',x,E-Ix)
box off
legend('boxoff')
legend('stimulus, E','response, R','location','northwest')
xlabel('x')



