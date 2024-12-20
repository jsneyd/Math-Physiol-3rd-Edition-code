%-------------------------------------------------------------------
%
% Matlab code for plotting the flux through a glucose transporter.
%
% For Chapter 2, Fig. 2.7 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
% Written by James Keener and James Sneyd
%
%-------------------------------------------------------------------

function glucose_transporter_plots
clear all
close all
clc

set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

k = 0.5;
 se=[0:.1:10];
si=0;
J0 = flux(se,si,k);
si=1;
J1 = flux(se,si,k);
si=2;
J2 = flux(se,si,k);
figure(1)
plot(se,J0,'g',se,J1,'r',se,J2,'b',se,zeros(length(se),1),'--','linewidth',2)
legend('\sigma_i=0','\sigma_i=1','\sigma_i=2','fontsize',18)
xlabel('External glucose','fontsize',20)
ylabel('Glucose flux','fontsize',20)

end % of main

function J=flux(se,si,k)
 J=(se-si)./((si+1+k).*(se+1+k)-k^2);

end
