
% Keener and Sneyd. Mathematical Physiology, 3rd edition, Chapter 14,
% Figure 1.4.

% The graph should be computed in concentrations and then converted to
% pressures using the solubility. Be careful with the units.

% You need to be careful with the ode solution otherwise you get wibbles
% and wobbles around the final steady state.

clearvars
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0);

global p

p.sigmaO = 1.4e-6;      % Units of Molar/mm Hg
p.KO = 30*p.sigmaO;
p.PO2 = 104*p.sigmaO;   % convert pressure to conc (units of Molar)
p.Z0 = 0.002;           % 2 mM

L = 100;
xspan = linspace(0,L,200);
init = 40*p.sigmaO;

opts = odeset('RelTol',1e-7,'AbsTol',1e-7,'MaxStep',0.01);
[X,W] = ode15s(@rhs,xspan,init,opts);
plot(X,W/p.sigmaO,'r')          % plot the pressure, not the conc
xlabel('dimensionless distance')
ylabel('[O_2]')


%% 
function out = rhs(x,W)
global p
fprime = 4*W^3*p.KO^4/( (p.KO^4 + W^4)^2 );
out = (p.PO2 - W)/((1+4*p.Z0*fprime));
end


