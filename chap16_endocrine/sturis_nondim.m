function main
clear all;

global L1 L2 L3 L4 A1 A2 A3 B0 B1

L1=0.4; L2=0.27; L3=0.06; L4=0.167; A1=21.7; A2=-39; A3=18; B0=-0.016; B1=0.06;

initial = [1 1 1 1 1 1];
tspan = linspace(0,200,1000);

[tt,sol]=ode45(@sturisfun,tspan,initial);

figure(1)
plot(tt,sol(:,1))
figure(2)
plot(tt,sol(:,3))  % Output per dL, as in Sturis
%save sturis.dat sol -ASCII

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt=sturisfun(t,y)
global L1 L2 L3 L4 A1 A2 A3 B0 B1

xx=y(1); yy=y(2); zz=y(3); h1=y(4); h2=y(5); h3=y(6); 

dydt = [input(t) - L1*(xx-L2*yy) - xx
        L1*(xx-L2*yy) - L3*yy
        quintic(h3) - B0*zz - B1*yy*zz
        L4*(xx - h1)
        L4*(h1-h2)
        L4*(h2-h3)
        ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ii=input(t)
ii=1 + 0.9*sin(t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=quintic(x)
global L1 L2 L3 L4 A1 A2 A3 B0 B1
out=1 + A1*x + A2*x*x + A3*x*x*x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

