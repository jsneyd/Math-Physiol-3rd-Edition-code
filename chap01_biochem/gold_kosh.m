%-------------------------------------------------------------------

% Matlab code for plotting Goldbeter-Koshland functions.

% For Chapter 1, Fig. 1.6 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

close all
clear all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);

w = linspace(0,1,100);

K1vals = [0.1,1.1,0.1];
K2vals = [0.05,1.2,1.2];
 
K1=0.1;
K2 = 0.05;
vrat = (1-w).*(K1+w)./((w.*(K2+1-w)));
plot(vrat,w,'r','LineWidth',2)
hold on
xlim([0,4])
 
K1=1.1;
K2 = 1.2;

vrat = (1-w).*(K1+w)./((w.*(K2+1-w)));
plot(vrat,w,'b--','LineWidth',2)
hold on
xlim([0,4])
 
K1=0.1;
K2 = 1.2;

vrat = (1-w).*(K1+w)./((w.*(K2+1-w)));
plot(vrat,w,'g--','LineWidth',2)
hold on
xlim([0,4])
 legend('boxoff')
  legend('K_1=0.1, K_2=0.05','K_1=1.1, K_2=1.2','K_1=0.1, K_2=1.2')
 box off
xlabel('v_1/v_2')
ylabel('w')