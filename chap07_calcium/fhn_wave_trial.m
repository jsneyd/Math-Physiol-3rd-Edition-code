%integrate the FHN traveling wave equations
function fhn_wave_analysis
close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 

%parameters

p.Iapp=0.01
p.alpha=0.1;
p.s=1;
p.eps=1;
p.gamma=0.5 ;

 init = [0,0,0];
dt=0.1;
tend=10;

tspan = [0:dt:tend];
[t,sol] = ode15s(@(t,x)rhs(t,x,p),tspan,init);

figure(1)
plot(t,sol(:,1),'r',t,sol(:,2),'b')  % time series
xlabel('Time')
legend('u','v')

% define the rhs here
function out= rhs(t,x,p)

u=x(1);  
v=x(2); %IPR inhibition
w = x(3);


 
out(1) = v;
out(2) = -p.s*v/p.eps-(u*(u-p.alpha)*(1-u)-w+p.Iapp)/p.eps^2;
out(3) = (u-p.gamma*w)/p.s

out = out';
 