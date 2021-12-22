function main
clear all;

global V1 t1 V2 t2 V3 t3 E

V1=3; t1=6; V2=11; t2=100; V3=10; t3=12; E=0.2;

initial = [1 1 1 1 1 1];
tspan = linspace(0,2000,1000);

[tt,sol]=ode45(@sturisfun,tspan,initial);

figure(1)
plot(tt,sol(:,1)/V1)
figure(2)
plot(tt,sol(:,3)/(10*V3))  % Output per dL, as in Sturis

save sturis.dat sol -ASCII

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt=sturisfun(t,y)
global V1 t1 V2 t2 V3 t3 E

xx=y(1); yy=y(2); zz=y(3); h1=y(4); h2=y(5); h3=y(6); 

dydt = [eff1(zz) - (xx/V1-yy/V2)*E - xx/t1
        (xx/V1-yy/V2)*E - yy/t2
        eff5(h3) + input(t) - eff2(zz) - eff3(zz)*eff4(yy)
        (xx - h1)/t3
        (h1-h2)/t3
        (h2-h3)/t3
        ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ii=input(t)
ii=216;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=eff1(z)
global V1 t1 V2 t2 V3 t3 E
out=209/(1+exp(-z/(300*V3)+6.6));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=eff2(z)
global V1 t1 V2 t2 V3 t3 E
out=72*(1-exp(-z/(144*V3)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=eff3(z)
global V1 t1 V2 t2 V3 t3 E
out=0.01*z/V3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=eff4(z)
global V1 t1 V2 t2 V3 t3 E
out=4 + 90/(1+exp(-1.772*log(z*(1/V2+1/(E*t2))) + 7.76));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=eff5(z)
global V1 t1 V2 t2 V3 t3 E
out=180/(1+exp(0.29*z/V1-7.5));