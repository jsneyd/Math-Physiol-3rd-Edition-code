
%  -------------------------------------------------------------------
%
%  Plot solutions of the lac operon model 
%
%   For Chapter 10, Section 10.2.5 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd
%
%  -------------------------------------------------------------------

close all
clear all
clc

set(0,                           ...
    'defaultaxesfontsize', 20,   ...
    'defaultaxeslinewidth', 2.0, ...
    'defaultlinelinewidth', 2.0);
global tmax
tmax = 10000;
nT = 1001;
init=[0.01 0.01 0.01 0.0 0];
tspan = linspace(0,tmax,nT );
 
 opts = odeset('AbsTol',1e-9,'RelTol',1e-9,'MaxStep',0.01); % be careful with the numerics
[T,Y] = ode15s(@(t,x)rhs(t,x),tspan,init,opts);

for i=1:nT 
    Le(i) = getLe(T(i));  % super inefficient, but it's fast enough and it works
end

figure(1)
    plot(T,Y(:,4),'LineWidth',2)
    ylabel('A (mM)')
    yyaxis right
    plot(T,Le,'LineWidth',2)
ylabel('L_e (mM)')
box off

figure(2)
semilogy(Le,Y(:,4))
xlabel('L_e')
ylabel('A')
box off


 %%
function out = rhs(t,x)
M = x(1);
B = x(2);
P = x(3);
A = x(4);
L = x(5);

aA=1.76e4; bA=2.15e4; gamA=0.52;
aB=1.66e-2; gamB=2.26e-2;
aP=10.0; gamP=0.65;
aM=9.97e-4; gamM=0.41;
aL=2880.0; gamL=2.26e-2;
Kaa=2.52e4; K=6000.0; KA=1.95; KL=9.7e-7; KLe=0.26; KL1=1.81;

Le = getLe(t);

out(1)=aM*(1.0+Kaa*(A)^2)/(K+Kaa*(A)^2) - (gamM)*M;
out(2)=aB*M-(gamB)*B;
out(3)=aP*M-(gamP)*P;
out(4)=aA*B*L/(KL+L)-bA*B*A/(KA+A)-(gamA)*A;
out(5)=aL*P*(Le)/(KLe+Le) - aA*B*L/(KL+L) - (gamL)*L;

out = out';
end

%% 
function Le = getLe(t)
global tmax
tmax

if t<=tmax/2
    Le = 0.04*t./(0.1*tmax+t);
else
    Le = 0.02*tmax/(0.1*tmax+tmax/2)./(1+(t-tmax/2)/(tmax/4));
    %Le=1./(1+(t-tmax/2)/(tmax/4));
end
end