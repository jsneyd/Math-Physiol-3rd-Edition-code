% Compute solution of the FHN equations.
% For Keener and Sneyd, Chapter 12

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 

p.eps = 0.02;  
p.alpha = -0.05; 
p.gamma = -0.6;
p.Iapp = 0;
dt = 0.001;
tend=7;

init = [0.2,0.49];
tspan = [0:dt:tend];
[t,sol] = ode15s(@(t,x)fhnrhs(t,x,p),tspan,init);

figure(1)
plot(t,sol(:,1),'r',t,sol(:,2),'b')  % time series
xlabel('Time'); ylabel('v');

% ----------------
function out=fhnrhs(t,x,p)

v=x(1);
w=x(2);
out(1) = (v.*(p.alpha-v).*(v-1) - w + p.Iapp)/p.eps;
out(2) = v-p.gamma*w;
out = out';
end