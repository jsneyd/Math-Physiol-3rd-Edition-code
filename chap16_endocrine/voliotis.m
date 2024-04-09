%    -------------------------------------------------------------------
%
%     Model of GnRH pulse generator by Voliotis et al.
%
%     For Chapter 16, Section 16.2.3 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function voliotis

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0);

par.I0 = 0.1;
par.dD = 0.25;
par.dN = 1;
par.dv = 10;
par.kkD = 4.5;
par.kkN = 320;
par.kD0 = 0.175;
par.kN0 = 0;
par.pv = 1;
par.v0 = 30000;
par.KD = 0.3;
par.KN = 32;
par.Kv1 = 1200;
par.Kv2 = 1200;

I0vals = [0,0.4/5,0.2,0.8];
for j= 1:length(I0vals)
    par.I0 = I0vals(j);
    tspan = linspace(0,140,1000);
    u0 = [ 0,0,0];
    [tout,U]=ode15s(@(t,u)derivs(t,u,par),tspan,u0);
    D = U(:,1);
    N = U(:,2);
    v = U(:,3);
    figure(1)
    plot(tout,v)
    hold on
end
hold off
axis([60 140 0 4000])
legend('0 Hz','2 Hz','5 Hz','20 Hz')
box off
xlabel('time (min)')
ylabel('v (firing rate)')

end % of main

%% RHS of ode
function out = derivs(t,u,par)

    D = u(1);
    N = u(2);
    v = u(3);

    fD = par.kD0 + par.kkD*v^2/(v^2+par.Kv1^2);
    fN = par.kN0 + par.kkN*(v^2/(v^2+par.Kv2^2))*(par.KD^2/(par.KD^2 + D^2));
    I = par.I0 + par.pv*v*N^2/(N^2 + par.KN^2);
    fv = par.v0*(1-exp(-I))/(1+exp(-I));

    Dp = fD - par.dD*D;
    Np = fN - par.dN*N;
    vp = fv - par.dv*v;

    out = [Dp;Np;vp];
end
