
%  -------------------------------------------------------------------
%
%   Ottesson model of oscillations in the baroreceptor loop.
%
%   For Chapter 11, Section 11.6.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------


clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global n tau ca R Vs alpha  beta mu

n = 7;
tau= 2;

ca=1.55;
R=1.05;
Vs=67.9 ;
alpha=0.84;
beta=1.17;
mu=93.0;

P0=88.71759;
H0=0.8026246;

%specify the output points
tspan = [0:.01:200];

%specify initial data
lags = tau;
init=[P0,H0]; %init is not used

sol = dde23(@ddeRHS,lags, @history, tspan);

figure(1)
plot(sol.x,sol.y(1,:),'r',sol.x,sol.y(2,:),'b','linewidth',2)
xlabel('t','fontsize',16)
legend('boxoff')

set(gca,'linewidth',1.5)
box off
legend('P','H')

figure(2)
plot(sol.y(2,:), sol.y(1,:))
xlabel('H')
ylabel('P')



%%
function out = ddeRHS(t,Y,Z)
    global   ca R Vs alpha  beta
    P=Y(1);
    H = Y(2);
    ylag = Z(1);

    dPdt = -P/(ca*R) + (Vs/ca)*H;
    dHdt = alpha*gs(ylag) - beta*(1-gs(P));
    out=[dPdt;dHdt];
    end

 %%
function gsout = gs(x)
    global mu n
    gsout = mu^n/(mu^n+x^n);
    end

%%
function s = history(t)
    s = ones(2,1);
    end

