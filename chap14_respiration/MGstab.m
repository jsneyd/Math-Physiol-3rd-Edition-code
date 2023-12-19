% MG stability plot

clearvars
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
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

