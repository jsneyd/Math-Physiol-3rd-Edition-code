% Compute solution of the FHN equations.
% For Keener and Sneyd, Chapter 5

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 14,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0); 

global p
p.eps = 0.01; p.alpha = 0.1; p.gamma = 0.5;
p.Iapp = 0.5;

tend = 4; nt = 50000;
init = [0.2,0];
tspan = linspace(0,tend,nt);
[t,sol] = ode15s(@(t,x)fhnrhs(t,x),tspan,init);

figure(1)
plot(t,sol(:,1),'r')  % time series


figure(2)  % phase plane
v = linspace(-0.4,1.4,100);
w1 = v.*(p.alpha-v).*(v-1) + p.Iapp;
w2 = v/p.gamma;
plot(sol(:,1),sol(:,2),'r',v,w1,'g--',v,w2,'b--')
ylim([-0.1,1.2*max(w1)]);
xlabel('v'); ylabel('w');


% ----------------
function out=fhnrhs(t,x)
global p
v=x(1);
w=x(2);
out(1) = (v.*(p.alpha-v).*(v-1) - w + p.Iapp)/p.eps;
out(2) = v-p.gamma*w;
out = out';
end