% phase portrait for blocking propagation failure
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)
R1 = 0.04;
R2 = 1.0;
a  = 0.25


u = [0:.01:1];

F = -u.^4/4 + ((a + 1)*u.^3)/3 - a*u.^2/2;
% curve 1:
 V1=sqrt(-F/R1);

 F1 = F(end);
V2=sqrt((F(end)-F)/R2);

plot(u,V1,u,V2,'--')
text(0.43,0.05,'i=1','fontsize',18)
text(0.75,0.15,'i=2','fontsize',18)
xlabel('V','fontsize',20)
ylabel('I','fontsize',20)
box off


