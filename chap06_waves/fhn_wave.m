
%  -------------------------------------------------------------------
%
%   Use the method of lines to compute the traveling wave in the
%   FitzHugh-Nagumo equations.
%
%   For Chapter 6, Exercise 6.7 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 


close all
clear all
clc
global n lam

set(0,                           ...
    'defaultaxesfontsize', 20,   ...
    'defaultaxeslinewidth', 2.0, ...
    'defaultlinelinewidth', 2.0);
n = 200;
L = 150;
diff = 1;
delx = L/(n-1);
lam = diff/(delx^2);
x = linspace(0,L,n);

v0 = zeros(1,n);   % initial condition
w0 = zeros(1,n);
v0(1:15) = 1;
%v0(end-15:end)= 1;    % for starting a wave at each end
tspan = linspace(0,200,10);
[t,sol] = ode15s(@(t,v)bs_mol(t,v),tspan,[v0 w0]);
tout = [2,5,8,10;];   % the output times
plot(x,sol(tout(1),1:n),x,sol(tout(2),1:n),x,sol(tout(3),1:n),x,sol(tout(4),1:n))
xlabel('x')
ylabel('v')
formatSpecF = '%3.0f\n';

legend(strcat('t = ',sprintf(formatSpecF,t(tout(1)))),strcat('t = ',sprintf(formatSpecF,t(tout(2)))),...
    strcat('t = ',sprintf(formatSpecF,t(tout(3)))),strcat('t = ',sprintf(formatSpecF,t(tout(4)))),'fontsize',16)


%%
function out = bs_mol(t,x)
global n lam

v = x(1:n);
w = x(n+1:2*n);
% no-flux boundary condition
out(1) = lam*(-2*v(1) +  2*v(2) ) + react(v(1),w(1));
out(2:n-1) = lam*(v(3:n) - 2*v(2:n-1) + v(1:n-2)) + react(v(2:n-1),w(2:n-1));
out(n) = lam*(-2*v(n) + 2*v(n-1)) + react(v(n),w(n));

out(n+1:2*n) = 0.01*(v(1:n) - w(1:n));   % The w equations

out = out';

end

%%
function out = react(v,w)
out = v.*(0.1-v).*(v-1) - w;
end