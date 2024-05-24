%
%  -------------------------------------------------------------------
%
%  This is an ode integrator for the FitzHugh-Nagumo equations.
%
%   For Chapter 5 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd
%
%  -------------------------------------------------------------------
function fhn

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

p.eps = 0.01;
p.alpha = 0.1;
p.gamma = 0.5;
p.Iapp = 0;

dt = 0.001;
tend=2;

init = [0.2,0];
tspan = [0:dt:tend];
[t,sol] = ode15s(@(t,x)fhnrhs(t,x,p),tspan,init);

figure(1)   % time series
plot(t,sol(:,1),'r',t,sol(:,2),'b')  % time series
xlabel('Time'); ylabel('v');
box off

figure(2)  % phase plane
v = linspace(-0.4,1.4,100);
w1 = v.*(p.alpha-v).*(v-1) + p.Iapp;
w2 = v/p.gamma;
plot(v,w1,'g--',v,w2,'b--',sol(:,1),sol(:,2),'r')
ylim([-0.1,1.2*max(w1)]);
xlabel('v'); ylabel('w');
box off
legend('boxoff')
legend('dv/dt=0','dv/dt=0')

end % of main


%% the differential equations

function out=fhnrhs(t,x,p)
v=x(1);
w=x(2);
out(1) = (v.*(p.alpha-v).*(v-1) - w + p.Iapp)/p.eps;
out(2) = v-p.gamma*w;
out = out';
end
