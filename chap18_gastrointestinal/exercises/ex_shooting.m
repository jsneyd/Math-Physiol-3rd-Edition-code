function stand_grad_shooting

clear all
close all
clc

L= 20;
cinit = 0.7;  %0.6169954  is the value needed to get c(L)=c0;
y0 = [cinit 0 0];
xspan = linspace(0,L,500);
[x,y] = ode15s(@(x,y)oderhs(x,y,L),xspan,y0);

ax = plotyy(x,y(:,1),x,y(:,3));
xlabel ("x",'FontSize',20);
ylabel (ax(1), "c",'FontSize',20);
ylabel (ax(2), "v",'FontSize',20);
end

%% 
function dydx = oderhs(x,y,L)
D = 1000;
r = 0.05;
c0 = 0.3;
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