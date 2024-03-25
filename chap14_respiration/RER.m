%  -------------------------------------------------------------------
%
%  Plot the respiratory exchange ratio.
%
%   For Chapter 14, Sections 14.5.1 and 14.5.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function RER

clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0);

global KO2

THb = 2.2e-3; %M

sigmaO2 = 1.4e-6;
sigmaCO2 = 3.3e-5;
PvO2 = 40.0;
PvCO2 = 45.0;
PiO2 = 150.0;
RT = 1.7e4;
Z0 = 2.2e-3;
KO2 = 30*sigmaO2;
Kc = 12.0;
Wv = (sigmaO2*PvO2);
phi = sigmaCO2*RT*(1.0+Kc);
phi1 = sigmaCO2*(1.0+Kc);

PaO2 = linspace(40,150,100);
VQO2 = (RT./(PiO2-PaO2)).*( sigmaO2*(PaO2-PvO2) + 4*Z0*( hemo(sigmaO2*PaO2) - hemo(Wv) ) );
PaCO2 = phi*PvCO2./(VQO2+phi);

% the RER curve
plot(PaO2,PaCO2,'LineWidth',2)
ylabel('P_{a,CO_2} (mm Hg)')
xlabel('P_{a,O_2} (mm Hg)')
hold on


% plot the separate blood and gas RER curves

Rgas = 0.8;
Pag1 = Rgas*(PiO2-PaO2);
plot(PaO2,Pag1,'r','LineWidth',1)
Rgas = 1.8;
Pag2=Rgas*(PiO2-PaO2);
plot(PaO2,Pag2,'--r','LineWidth',1)

Rblood=0.8;
Pab1 = (sigmaCO2*(1+Kc)*PvCO2 -  Rblood*( sigmaO2*(PaO2-PvO2) + 4*Z0*( hemo(sigmaO2*PaO2) - hemo(Wv) ) ) )/(sigmaCO2*(1+Kc));
plot(PaO2,Pab1,'g','LineWidth',1)
Rblood=1.8;
Pab2 = (phi1*PvCO2 -  Rblood*( sigmaO2*(PaO2-PvO2) + 4*Z0*( hemo(sigmaO2*PaO2) - hemo(Wv) ) ) )./phi1;
plot(PaO2,Pab2,'--g','LineWidth',1)

ylim([0,50])

end


%%

function out = hemo(PO)
global KO2
out = PO.^4./(KO2^4+PO.^4);
end