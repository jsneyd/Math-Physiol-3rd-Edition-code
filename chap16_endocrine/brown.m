function brown
close all
clear all
clc

initial = [0.1 0 0];
tspan = linspace(0,15,100000);
[tt,sol]=ode15s(@brownfun,tspan,initial);

figure(1)
plot(tt,sol(:,1))
figure(2)
plot(tt,sol(:,3))

save brown.dat sol -ASCII

%%
function dydt=brownfun(t,y)
s0=50; c=0.2; k1=1; k2=0.02; k3=0.02; k4=1; rg=2.5;
p1=100; p2=100; p3=0.3; rz=1;

v=y(1); g=y(2); z=y(3);
h = max([v 0]);
p = p1/(1+exp(-p2*(h-p3)));
a = getstim(t);
dydt = [s0*(-v*(v-c)*(v-1) - k1*g + k2*a)
        k3*a + k4*v - rg*g
        p - rz*z
        ];
%%
function a=getstim(t)
a=0;
stimperiod = 1;
stimwidth = 0.9*stimperiod;
if ( rem(t,stimperiod)<stimwidth) 
    a=4;
end