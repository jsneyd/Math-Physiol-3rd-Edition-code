% Compute solution of the FHN equations.
% For Keener and Sneyd, Chapter 5

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 


p.eps = 0.01; 
p.alpha = 0.1; 
p.gamma = 0.5;
% pick a value of p.Iapp
p.Iapp = 0.5;
p.Iapp = 0.;

dt = 0.01;
if(p.Iapp == 0)
    tend=2;
else
    tend=4;
end

 
init = [0.2,0];
tspan = [0:dt:tend];
[t,sol] = ode15s(@(t,x)fhnrhs(t,x,p),tspan,init);

figure(1)
plot(t,sol(:,1),'r')  % time series
xlabel('Time'); ylabel('v');

figure(2)  % phase plane
v = linspace(-0.4,1.4,100);
w1 = v.*(p.alpha-v).*(v-1) + p.Iapp;
w2 = v/p.gamma;
plot(v,w1,'g--',v,w2,'b--',sol(:,1),sol(:,2),'r')
ylim([-0.1,1.2*max(w1)]);
xlabel('v'); ylabel('w');
legend('dv/dt=0','dn/dt=0')

% ----------------
function out=fhnrhs(t,x,p)

v=x(1);
w=x(2);
out(1) = (v.*(p.alpha-v).*(v-1) - w + p.Iapp)/p.eps;
out(2) = v-p.gamma*w;
out = out';
end