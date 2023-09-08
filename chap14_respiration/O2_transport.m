
% Keener and Sneyd. Mathematical Physiology, 3rd edition, Chapter 14,
% Figure 1.4.

% The graph should be computed in concentrations and then converted to
% pressures using the solubility. Be careful with the units. You really do
% not want to work in units of Molar. The values for sigma are way too
% small, and the numerics starts going a little weird. So I converted all
% units to microMolar. The code plays much more nicely.


clearvars
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0);

global p

%% The figure in the text

p.sigmaO = 1.4;      % Units of uMolar/mm Hg
p.KO = 30*p.sigmaO;
p.PO2 = 104*p.sigmaO;   % convert pressure to conc (units of uMolar)
p.THb = 2000;           % 2 mM

L = 100;
xspan = linspace(0,L,2000);
init = 40*p.sigmaO;

[X,O] = ode15s(@rhs,xspan,init);
plot(X,O/p.sigmaO,'r')          % plot the pressure, not the conc
xlabel('dimensionless distance')
ylabel('[O_2]')
hold on

%% This next bit is for one of the Exercises.

% Now do the solution without the equilibrium approximation. When solving
% this way you need to be careful with your initial conditions, as you need
% to start on the critical manifold. Otherwise, the solution takes some
% time to catch up and you get an offset in the solution. 

timescale = 10;
p.k1 = timescale;
p.km1 = p.KO^4*p.k1;
Oinit = 40*p.sigmaO;
Yinit = p.THb*Oinit^4/(p.KO^4 + Oinit^4);
init = [Oinit, Yinit];
opts = odeset('RelTol',1e-7,'AbsTol',1e-7,'MaxStep',0.01);
[X,sol] = ode15s(@twodrhs,xspan,init,opts);
O = sol(:,1); Y = sol(:,2);
plot(X,O/p.sigmaO,'--b')          % plot the pressure, not the conc
xlabel('dimensionless distance')
ylabel('[O_2]')


%% 
function out = rhs(x,O)
global p
fprime = 4*O^3*p.KO^4/( (p.KO^4 + O^4)^2 );
out = (p.PO2 - O)/((1+4*p.THb*fprime));
end

%%
function out = twodrhs(x,dum)
global p
O = dum(1); Y = dum(2);
X = p.THb - Y;
out(1) = p.PO2 - O + 4*p.km1*Y - 4*p.k1*X.*O.^4;
out(2) = p.k1*X*O.^4 - p.km1*Y;
out = out';
end
