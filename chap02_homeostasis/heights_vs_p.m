%-------------------------------------------------------------------

% Matlab code for plotting the heights of two water columns.

% For Chapter 2, Fig. 2.15 of
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


h1=[0:.01:2];
p=h1.^2-h1;
h2=2-h1;

figure(1)
plot(p,h1,p,h2)
xlabel('p')
legend('boxoff')
legend('\eta_1','\eta_2','location','northwest')
axis([0 2 0 2])
set(gca,'linewidth',1.5)
box off