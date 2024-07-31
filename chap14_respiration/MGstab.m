%  -------------------------------------------------------------------
%
%  Mackey-Glass stability plot.
%
%   For Chapter 14, Section 14.6 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

clearvars
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
y = [0:.01:3];
n=3;

F = y.^n./(1+y.^n);
Fp = n*y.^(n-1)./((1+y.^n).^2);
g = F-y.*Fp;

figure(1)
plot(y,F,y,Fp,y,g,y,zeros(1,length(y)),'--')
xlabel('y')
text(1.4,.86,'F(y)','fontsize',18)
text(2,.2,'dF/dy','fontsize',18)
text(1.7,.6,'g(y)','fontsize',18)

box off

