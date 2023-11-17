%-------------------------------------------------------------------

% Matlab code for plotting the velocity of a cooperative enzyme.

% For Chapter 1, Fig. 1.4 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

close all
clear all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

e0 = 1; k2 = 1; k4 = 2;

s = linspace(0,2,100);
K1 = 1000; K2 = 0.001;
V = (k2*K2 + k4*s)*e0.*s./(K1*K2 + K2*s + s.^2);
plot(s,V,'r','LineWidth',2)
hold on

K1 = 0.5; K2 = 2;
V = (k2*K2 + k4*s)*e0.*s./(K1*K2 + K2*s + s.^2);
plot(s,V,'--b','LineWidth',2)

K1 = 0.5; K2 = 100;
V = (k2*K2 + k4*s)*e0.*s./(K1*K2 + K2*s + s.^2);
plot(s,V,'--g','LineWidth',2)

xlabel('s')
ylabel('V')