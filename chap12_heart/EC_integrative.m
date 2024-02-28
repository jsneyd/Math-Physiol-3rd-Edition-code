%  -------------------------------------------------------------------
%
%   A simple integrative model, combining multiple CRUs to get a whole-cell
%   response.
%
%   For Chapter 12, Section 12.3.3 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function EC_integrative

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);


% The variables are, in order:
% op - open probability of the L-type channel
% n3 - probability the L-type channel is in state N
% cd - Ca conc in the diadic cleft
% Obar - open probability of the RyR
% c - Ca conc in the myoplasm.


VGCC_mean = 0.05;               % mean L-type density
VGCC_sig = 0.04;                % spread of L-type densities
kRyR_mean = 20;                 % RyR density
kRyR_sig = 15;                  % spread of RyR densities

% RyR parameters
p.K1 = 0.55;
p.k2 = 0.05;
p.km2 = 0.01;

% calcium parameters
p.Vc = 0.2;          % PM pump density
p.Kc = 0.2;
p.D = 10;
p.Vf = 0.01;

% L-type parameters
p.f=0.85;
p.omega=0.02;
p.g=2;


% Select VGCC and kRyR from a
% uniform distribution and then calculate the responses for all Vstim
Nsample = 100;
for j = 1:Nsample
    j                                       % just to track progress
    % first find the unstimulated steady state (at Vm = -75)
    p.Vm = -75;
    p.kL = 2*VGCC_sig*rand + VGCC_mean - VGCC_sig;
    p.kRyR = 2*kRyR_sig*rand + kRyR_mean - kRyR_sig;
    dt = 0.01;
    tend=250;
    init = [0.0041, 0.4, 0.006469, 0.75559, 0.006];
    tspan = [0:dt:tend];
    [t,sol] = ode15s(@(t,x)Ltrhs(t,x,p),tspan,init);
    init = sol(end,:);                      % the steady state for all the following Vstim runs

    % Now start at the steady state and stimulate by increasing Vm
    n = 150;                                % number of voltage-clamp stimulations
    Vstim = linspace(-74,65,n);             % stimulating voltage clamp levels

    for i = 1:n
        p.Vm = Vstim(i);
        dt = 0.01;
        tend=10;
        tspan = [0:dt:tend];
        [t,sol] = ode15s(@(t,x)Ltrhs(t,x,p),tspan,init);
        op = sol(:,1);
        cd = sol(:,3);
        Obar = sol(:,4);

        % pull out the peak fluxes for plotting
        peak_JRyR(j,i) = max(RyR_flux(Obar,cd,p));
        peak_Ltype(j,i) = max(-Ltype(op,p));
    end
end

peak_JRyR_mean = mean(peak_JRyR,1);
peak_Ltype_mean = mean(peak_Ltype,1);
plot(Vstim, peak_JRyR_mean, Vstim, peak_Ltype_mean)
hold on


% Finally, plot the Vstim curve for the parameters in the middle of the
% distribution. Easiest, but less efficient, to recalculate it.

% first find the unstimulated steady state (at Vm = -75)
p.Vm = -75;
p.kL = VGCC_mean;
p.kRyR = kRyR_mean;
dt = 0.01;
tend=250;
init = [0.0041, 0.4, 0.006469, 0.75559, 0.006];
tspan = [0:dt:tend];
[t,sol] = ode15s(@(t,x)Ltrhs(t,x,p),tspan,init);
init = sol(end,:);                      % the steady state for all the following Vstim runs

% Now start at the steady state and stimulate by increasing Vm
n = 150;                                % number of voltage-clamp stimulations
Vstim = linspace(-74,65,n);             % stimulating voltage clamp levels

for i = 1:n
    p.Vm = Vstim(i);
    dt = 0.01;
    tend=10;
    tspan = [0:dt:tend];
    [t,sol] = ode15s(@(t,x)Ltrhs(t,x,p),tspan,init);
    op = sol(:,1);
    cd = sol(:,3);
    Obar = sol(:,4);

    % pull out the peak fluxes for plotting
    peak_JRyR_fix(i) = max(RyR_flux(Obar,cd,p));
    peak_Ltype_fix(i) = max(-Ltype(op,p));
end
plot(Vstim, peak_JRyR_fix, Vstim, peak_Ltype_fix,'LineWidth',2)

%save EC_integrated.mat Vstim peak_JRyR_mean peak_JRyR_fix peak_Ltype_mean peak_Ltype_fix

end % of main

%%
function out=Ltrhs(t,x,p)

    op=x(1);
    n3=x(2);
    c3=1-n3-op;
    cd = x(3);
    Obar = x(4);
    c = x(5);

    p.alpha = 2*exp(0.012*(p.Vm-35));
    p.beta = 0.0882*exp(-0.05*(p.Vm-35));
    p.gamma = 0.44*cd;
    p.fbar = p.alpha*p.f/(p.alpha+p.beta);
    p.omegabar = p.omega*(p.beta/2+p.alpha)/(2*p.alpha+p.beta/2);
    p.gammabar = p.gamma*(p.beta+2*p.alpha)/(p.alpha+p.beta);

    % Ca fluxes
    JLtype = Ltype(op,p);
    Jpm = p.Vc*c^2/(p.Kc^2 + c^2);
    JRyR = RyR_flux(Obar,cd,p);

    % L-type channel
    dopdt = -p.g*op+p.fbar*n3;
    dn3dt =-(p.fbar+p.gammabar)*n3 +p.omegabar*c3+p.g*op;

    % cd equation
    dcddt = -JLtype + JRyR - p.D*(cd - c);

    % RyR equation
    dObardt = p.km2*(1-Obar) - p.k2*cd*Obar;

    % cytosolic calcium equation
    dcdt = p.Vf*p.D*(cd-c) - Jpm;

    out=[dopdt; dn3dt; dcddt; dObardt; dcdt];
end

%% L-type current
function out = Ltype(op,p)
    out = p.kL*op*(p.Vm - 60); % Nernst potential is around 50 - 70 mV.
end

%% RyR flux
function out = RyR_flux(Obar,cd,p)
    out = p.kRyR*Obar.*(cd.^2)./(p.K1^2 + cd.^2);
end





