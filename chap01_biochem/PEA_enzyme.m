%-------------------------------------------------------------------

% Matlab code for simulating the equilibrium approximation (PEA)
% of enzyme kinetics.

% For Chapter 1 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

function enzyme
% asymptotic analysis of enzyme kinetics equations - PEA
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
 

% Part 1:  Equilibrium approximation
par.a = 2.2;
par.b = 1.7;
par.e = 0.1;

tspan = linspace(0,20,2000);
initial = [1,1];
[t,Y] = ode15s(@(t,y)rhseq(t,y,par),tspan,initial);
s=Y(:,1);
z=Y(:,2);
x = (z-s)/par.a;   % the original variable, x
figure(1)
plot(t,s,t,z,'linewidth',2)
legend('\sigma','z','fontsize',18)
xlabel('\tau','fontsize',20)
axis([0 2 0 1])

% now plot the slow manifold and the solution together
slow_s = linspace(0,1,100);
slow_z = slow_s + par.a*par.b*slow_s./(1+par.b*slow_s);
slow_x = (slow_z-slow_s)/par.a;

figure(2)
plot(s,z,slow_s,slow_z,'--','linewidth',2)
xlabel('\sigma','fontsize',20)
ylabel('z','fontsize',20)
figure(3)
plot(s,x,slow_s,slow_x,'--','linewidth',2)
xlabel('\sigma','fontsize',20)
ylabel('x','fontsize',20)

end

%%
function out = rhseq(t,y,par)
s = y(1);
z = y(2);
out(1) = z-s - par.b*s*(par.a + s - z);
out(2) = par.e*(s-z);
out = out';  % don't forget this or you get rubbish.
end
 