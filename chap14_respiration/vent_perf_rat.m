
%  -------------------------------------------------------------------
%
%  Plot the ventilation-perfusion ratio
%
%   For Chapter 14, Section 14.5 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function vent_perf_rat

clearvars
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

Pvco2=45;
Kc = 12;
RT = 1.7e4;
sigmaco2=3.3e-5; % mM/mm Hg
sigmao2=1.4e-6;
THb = 2.2e-3; %M


Paco2 = [0.1:.01:1]*Pvco2;

VbyQ=sigmaco2*RT*(1+Kc)*(Pvco2-Paco2)./Paco2;

PiO2 = 150;
PvO2 = 40;
PaO2 = [PvO2:.1:150];
Ta = sigmao2*PaO2+4*THb*f(PaO2);
Tv = sigmao2*PvO2+4*THb*f(PvO2);
VbyQO= RT*(Ta-Tv)./(PiO2-PaO2);

figure(1)
    plot(VbyQ,Paco2,VbyQO,PaO2)
    axis([0 2 30 130])
    xlabel('Ventilation-perfusion ratio')
    ylabel('Alveolar partial pressure (mm Hg)')
    text(1.5,43,'CO_2','fontsize',18)
    text(1.5,120,'O_2','fontsize',18)
    box off

end  % of main

function out = f(PO)
Ko2 = 30;
out = PO.^4./(Ko2^4+PO.^4);
end
