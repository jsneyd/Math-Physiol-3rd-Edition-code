%  -------------------------------------------------------------------
%
%   Calculate the steady state curves for the Positive feedback motif.
%
%   For Chapter 10, Section 10.1.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);

%Parameters
m = 8;
n = 4;
 
k1vals=[0.5,0.2,1.35];
xlim=[2.5,4,2];
Slim=[0.6,0.3,2];
k2 = 1;
for j = 1:3
    k1=k1vals(j);

x=[0:.01:xlim(j)];
% find y as a function of x;

f = x.^m./(1+x.^m);

y=f/k2;  

g = y.^n./(1+y.^n);
S=k1*x-g;

figure(j)
plot(S,x)
formatSpecF = '%6.2f\n';
 
 title(strcat('k_1 = ',sprintf(formatSpecF,k1)),'fontsize',18)
xlabel('S')
ylabel('x')
axis([0 Slim(j) 0 xlim(j)])

end
