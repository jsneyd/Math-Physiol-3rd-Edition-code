
%    -------------------------------------------------------------------
%
%     Hamer model of light adaptation in rods.
%
%     For Chapter 19, Section 19.3 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function hamer
clear all
close all
clc

global ng Jdark

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);

IC = [0,0,1.1068,0.7106,19.3133];
tspan = linspace(0,30,1000);
options = odeset('MaxStep',0.01);
ng = 3; Jdark = 72;

% First find the dark state, to use as initial condition
stim = 0;
[t,Y] = ode15s(@(t,X)rods(t,X,stim,0,0),tspan,IC,options);
IC = Y(end,:);


% compute flash responses
keep = tspan'; % for external plotting
ton = 1;
toff = 1.01;
stimlist = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000];
for i = 1:10
    stim = stimlist(i);
    [t,Y] = ode15s(@(t,X)rods(t,X,stim,ton,toff),tspan,IC,options);
    photocurrent = (Y(:,3)/IC(3)).^ng*Jdark;
    figure(1)
        plot(t,photocurrent(1) - photocurrent)
        hold on
    keep = [keep photocurrent(1) - photocurrent];  % for external plotting
end
xlabel('t (s)')
ylabel('Photocurrent (pA)')
box off

%writematrix(keep,'hamer.dat');   % for external plotting


% compute step responses
tspan = linspace(0,100,1000);
ton = 1;
toff = 60;
stimlist = [1e-8 2e-8 1e-7 2e-7 1e-6 1e-5 1e-4];
for i = 1:7
    stim = stimlist(i)*1e2;
    [t,Y] = ode15s(@(t,X)rods(t,X,stim,ton,toff),tspan,IC,options);
    photocurrent = (Y(:,3)/IC(3)).^ng*Jdark;
    figure(2)
        plot(t,photocurrent(1) - photocurrent)
        hold on
end
xlabel('t (s)')
ylabel('Photocurrent (pA)')
box off
end % of main



%% rhs for ODE

function out = rods(t,X,stim,ton,toff)
    R = X(1); E = X(2); g = X(3); c = X(4); b = X(5);
    global ng Jdark

    tauR = 0.4; tauE = 2;
    v = 1160; bt = 660;
    gamma = 73; kon = 0.1;
    koff = 0.8; cmin = 0.005;
    Kc = 0.8; nc = 3;
    Amax = 10; betaE = 1.7e2;
    beta_dark = 0.13;
    alpha = 0.8; gdark = 2;

    light = stim*(heaviside(t-ton) - heaviside(t-toff));
    f = (g/gdark)^ng;

    out(1) = light - R/tauR;
    out(2) = v*R - E/tauE;
    out(3) = Amax/(1+(c/Kc)^nc) - (beta_dark + betaE*E)*g;
    out(4) = alpha * f * Jdark - gamma*(c - cmin) - kon*(bt - b)*c + koff*b;
    out(5) = kon*(bt-b)*c - koff*b;

    out = out';
end
