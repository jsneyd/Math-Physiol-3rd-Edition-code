% This code determines the initial condition that gives the correct value
% for c at the end of the tube (in the standing-gradient model).

function ex_shooting

clear all
close all
clc

L= 100;  % cinit = 0.6169954  is the value needed to get c(L)=c0, if L=20.;
c0 = 0.3;
xspan = linspace(0,L,10000);

% Do a bisection method on cinit to find the best value

cup = 2;
cdown = 0.5;

for i=1:50
cinit = (cup + cdown)/2;
y0 = [cinit 0 0];
stop_cond = odeset('Events',@stopping);   % The stopping conditions for the integration
[x,y] = ode15s(@(x,y)oderhs(x,y,L,c0),xspan,y0,stop_cond);

% Now look at the end of the integration and thus decide how to move the
% endpoints of the bisection method.
if (y(end,1) > c0)
    cup = cinit;
else
    cdown = cinit;
end
plot(x,y(:,1))
ylim([0,2])
hold on
end

cinit    % plot out the final value of cinit, which is the value that gives c(L) = c0.

% ax = plotyy(x,y(:,1),x,y(:,3));
% xlabel ("x",'FontSize',20);
% ylabel (ax(1), "c",'FontSize',20);
% ylabel (ax(2), "v",'FontSize',20);

end

%% Define the ODEs
function dydx = oderhs(x,y,L,c0)
D = 1000;
r = 0.05;
P = 0.2;
N0 = 0.3;
a = 0.1;

N = N0*(1-heaviside(x-a*L));
c = y(1);
z = y(2);
v = y(3);

dydx(1) = z;
dydx(2) = (1/D)*(  2*P*c*(c-c0)/r + z*v - 2*N/r    );
dydx(3) = 2*P*(c-c0)/r;

dydx = dydx';
end

%% Define the condition under which we stop the integration because it's already going wrong
function [value, isterminal, direction] = stopping(x,y)
value = [y(1);y(1)-10];
isterminal = [1;1];
direction = [0;0];
end