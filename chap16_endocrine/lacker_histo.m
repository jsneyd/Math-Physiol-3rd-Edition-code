
function lacker_histo
clear all
close all
clc

for i=1:1500
    [t_ov(i),ov_num(i)] = lacker;
end
figure(1)
hist(t_ov,50)

figure(2)
hist(ov_num,3)

dlmwrite('lacker.dat',[t_ov' ov_num'])

end

%%

function [t_ov,ov_num] = lacker
par.n = 500;
par.M1 = 3;
par.M2 = 15;

fol0 = sort(0.1*rand(1,par.n),'descend'); % They have to be ordered.
xi10 = fol0(1);
gam0 = fol0/fol0(1);  % initial conditions for the gammas.
tspan = linspace(0,5,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% First solve for the gammas, xi1 and t. Easier to do it all together.
[tau,y]=ode15s(@(tau,y)rhs(tau,y,par),tspan,[gam0,xi10,0],options);
gam = y(:,1:par.n);
xi1 = y(:,par.n+1);
t = y(:,par.n+2);

% Now reconstruct and plot the xi
xi = gam.*xi1;
ov_index = find(xi(:,1)>2,1);

% Finally, pull out the ovulation time and number
t_ov = t(ov_index);
ov_num = sum(gam(end,:)>0.9);

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
