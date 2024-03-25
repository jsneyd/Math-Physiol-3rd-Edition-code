
%  -------------------------------------------------------------------
%
%   Code to simulate post-inhibitory rebound coupled HH neurons.
%
%   For Chapter 14, Section 14.7.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function PIR_HH
close all
clear all
clc
set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0);

global gpir   Vpir gL gsyn Vsyn phi thetasyn ksyn VL

gsyn = 0.3;
thetasyn = -44;
gpir = 0.3;
ksyn = 2;
Vpir = 120;
VL = -60;
Vsyn = -80;
gL = 0.1;
phi = 3;

% integrate the ode
init = [-50;1;-40;0]; %initial data for the odes
tstep = 0.1;
t_end = 300;

%specify the output points
tspan = [-50:tstep:t_end];

[T,S] = ode15s(@(t,x)rhs(t,x),tspan,init);
V1=S(:,1);
V2 = S(:,3);

figure(1)
    plot(T,V1,T,V2,'--','linewidth',2)
    legend('boxoff')
    legend('V_1','V_2')
    xlabel('Time')
    axis([0 300 -80 0])

end % of main

%% the right hand side for ode simulation: 

function out=rhs(t,s)
    global gpir  Vpir gL gsyn Vsyn phi VL
    
    V1 = s(1);
    h1 = s(2);
    V2 = s(3);
    h2 = s(4);
    
    out1 = gates(V1); 
    minf1 = out1(1);
    hinf1 = out1(2);
    tauh1 = out1(3);
    Sinf1 = out1(4);
    
    out2 = gates(V2);
    minf2 = out2(1);
    hinf2 = out2(2);
    tauh2 = out2(3);
    Sinf2 = out2(4);
    
    V1p = -gpir*(minf1^3)*h1*(V1-Vpir) - gL*(V1-VL) - gsyn*Sinf2*(V1-Vsyn);
    h1p = phi*(hinf1-h1)/tauh1 ;
    V2p = -gpir*(minf2^3)*h2*(V2-Vpir) - gL*(V2-VL) - gsyn*Sinf1*(V2-Vsyn);
    h2p = phi*(hinf2-h2)/tauh2 ;
    
    out = [V1p;h1p;V2p;h2p];
end

%%
function out = gates(V)
    global thetasyn ksyn
    minf = 1/(1+exp(-(V+65)/7.8));
    hinf = 1/(1+exp((V+81)/11));
    tauh=hinf*exp((V+162.3)/17.8);
    Sinf = 1/(1+exp(-(V-thetasyn)/ksyn));
    
    out = [minf,hinf,tauh,Sinf];
end