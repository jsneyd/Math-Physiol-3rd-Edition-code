%  -------------------------------------------------------------------
%
%  The compartmental model for oxygen uptake by the lungs.
%
%   For Chapter 14, Section 14.4.5 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function compartmental_O

clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0,...
   'defaultlinemarkersize',12);

global p

p.sigma = 1.4;
p.KO = 30*p.sigma;
p.THb = 2200;
p.RT = 1.7e-2;   % units of uM
p.D =  32.5;
p.PO =  150;  % air partial pressure
p.gamma = 37;

% solve the ODE for TO
init = 6000;
tspan = linspace(0,10,100);
[T,Y] = ode15s(@Tode,tspan,init);
for j = 1:length(T)
    [O(j),Otilde(j)] = getO_Otilde(Y(j));
end

figure(1)
    plot(T,Y)
    xlabel('time (nondimensional)')
    ylabel('[O_2] in Body')

% the target is Pa = 104, Pv = 40
Oa = 104*p.sigma;
conc(Oa)

figure(2)
    plot(T, O/p.sigma,'r', T , Otilde/p.sigma,'b')
    title('O_2 Partial Pressures')
    xlabel('time (dimensionless')
    O(end)/p.sigma
    Otilde(end)/p.sigma

end % of main



%% ODE for TO
function out = Tode(t,x)
    TO = x;
    global p
    [O,Otilde] = getO_Otilde(TO);
    out = p.D*(p.sigma*p.PO - Otilde) - p.gamma*O;
end

%% get O and Otilde as functions of TO
function [O,Otilde] = getO_Otilde(TO)
    global p
    options = optimset('Display','off');
    O = fsolve(@(x) conc(x)-TO,20,options);
    Otilde = fsolve(@(x) p.D*(p.sigma*p.PO - x) + TO - conc(x),20,options);
end

%% oxygen saturation function
function out=conc(O)
    global p
    out = O + 4*p.THb*O.^4./(p.KO^4 + O.^4);
end
