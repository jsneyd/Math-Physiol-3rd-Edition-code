
% Code to solve a simple version of the Longtin-Milton model of the pupil
% light reflex. This doesn't run under Octave, which hasn't yet implemented
% the delay differential equation solver, dde23.

function longtin_milton
clear all
close all
clc

options = ddeset('Reltol',0.00005);  % To get output at a decent resolution

lags = 1;
tspan = (0:10);
sol = dde23(@ddefun, lags, @history, tspan,options);

plot(sol.x,areafun(sol.y),'-o')
xlabel('Time (s)');
ylabel('Pupil area (mm^2)');
end

%%
function dydt = ddefun(t,y,Z)
gamma = 5;
I = 10;
phibar = 1;
ylag1 = Z(:,1);
dydt = gamma*capF(log(I*areafun(ylag1)/phibar))-y(1);
end

%%
function s = history(t)
  s = ones(1,1);
end

%% 
function out=capF(x)
out = heaviside(x)*x;
out = x.*(exp(x)).^1./(1+exp(x).^1);   % smooth approximation to the function with the kink. Many other possibilities.
end

%%
function out = areafun(x)
lam = 30;
n = 7;   % vary n to get different oscillatory behaviour
theta = 10;
out = lam*theta^n./(x.^n + theta^n);
end
