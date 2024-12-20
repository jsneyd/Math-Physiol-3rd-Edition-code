
% --------------------------------
% Code to plot the stability curve for  the Longtin-Milton model of the pupil
% light reflex.  

% Used to generate the image in Fig. 26 of Chapter 19 
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd
% --------------------------------

clear all
close all
clc
set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0, ...
'defaultpatchlinewidth', 0.7);

eta = [pi/2+.01:.01:pi];
G = -1./cos(eta);
delay = -eta./tan(eta);

figure(1)
    plot(delay,G)
    axis([0 10 0 10])
    xlabel('Delay')
    ylabel('Gain')
    box off
    text(1,1,'Stable','fontsize',18)
    text(4,4,'Unstable','fontsize',18)

    
    
