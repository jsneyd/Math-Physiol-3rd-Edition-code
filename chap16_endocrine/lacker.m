
%    -------------------------------------------------------------------
%
% Solve the Lacker model for the control of ovulation number. The model is
% solved in the transformed coordinates (gamma rather than xi) because it is
% much easier that way, as this avoids the problem of  solutions going to infinity in
% finite time, which is a real pain
%
%     For Chapter 16, Section 16.3.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function lacker

clear all; close all; clc;

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0);
par.n = 10;
par.M1 = 4;
par.M2 = 20;

fol0 = sort(0.1*rand(1,par.n),'descend'); % They must be ordered.
xi10 = fol0(1);
gam0 = fol0/fol0(1);  % initial conditions for the gammas.
tspan = linspace(0,2,1000);
options = odeset('RelTol',1e-7,'AbsTol',1e-7);

% First solve for the gammas, xi1 and t. Easier to do it all together.
[tau,y]=ode15s(@(tau,y)rhs(tau,y,par),tspan,[gam0,xi10,0],options);
gam = y(:,1:par.n);
xi1 = y(:,par.n+1);
t = y(:,par.n+2);
figure(1)
plot(tau,gam)
xlabel('\tau')
ylabel('\gamma_i')

% Now reconstruct and plot the xi
xi = gam.*xi1;
ov_index = find(xi(:,1)>5,1);
t_ov = t(ov_index)
ov_num = sum(gam(end,:)>0.9)

figure(2)
plot(t,xi)
ylim([0 5])
xlabel('t')
ylabel('\xi_i')


end
%%
function out=rhs(tau,y,par)
gam = y(1:par.n);
xi1 = y(par.n+1);
GG = sum(gam);
gam_deriv = gam.*(1-gam).*(par.M1*par.M2*(1+gam) - GG*(par.M1+par.M2));
xi_deriv = (1/xi1)*(1 - xi1^2*(GG-par.M1)*(GG-par.M2));
time_deriv = 1/(xi1^2);
out = [gam_deriv',xi_deriv,time_deriv]';
end
