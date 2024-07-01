%  -------------------------------------------------------------------
%
%   The trp operon.
%
%   For Chapter 10, Section 10.2.4 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);

krat=500;
konrat=5;
krrat=50;
ktrat=100;
mbyk1=1;
mbk2=2;

T=[0:.1:50];
R = T.^2./(ktrat+T.^2);

F = krat./(1+konrat+krrat*R);
 
figure(1)
plot(T,F,T,mbyk1*T,'--',T,mbk2*T,'--')

text(4,60,'F(T)','fontsize',20)
text(29,75,'\mu/K=2','fontsize',20)
text(40,36,'\mu/K=1','fontsize',20)
xlabel('T')

figure(2)
mbyK = F./T;
plot(mbyK,T)
ylabel('T')
xlabel('\mu/K')
axis([0 100 0 10])

