function main
clear all
close all
clc

initial = [0.01 0 0];
tspan = linspace(0,20,4000);

[tt,sol]=ode45(@brownfun,tspan,initial);

figure(1)
plot(tt,sol(:,1))
figure(2)
plot(tt,sol(:,3))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt=brownfun(t,y)
s0=200; c=0.2; k1=1; k2=0.02; k3=0.1; k4=1; k5=2.5;
b1=2; b2=0.59; b3=10; b4=0.1; 
p1=100; p2=100; p3=0.3; d1=1;

v=y(1); g=y(2); z=y(3);
b = b1 - b2/(1+exp(-b3*(v-b4)));
h = max([v v]);
p = p1/(1+exp(-p2*(h-p3)));

a = getstim(t);

dydt = [s0*(-v*(v-c)*(v-1) - k1*g + k2*a)
        b*(k3*a + k4*v - k5*g)
        p - d1*z
        ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a=getstim(t)

a=0;
stimperiod=5;
stimwidth=3;
if ( rem(t,stimperiod)<stimwidth) 
    a=3.5;
end