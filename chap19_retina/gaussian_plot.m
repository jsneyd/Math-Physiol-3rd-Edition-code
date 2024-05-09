% Matlab/Octave code to plot excitatory and inhibitory surround function
 
% For Figure 23, Chapter 19 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------

function gaussian_plot

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);
% parameters

s1 = 3;
s2 = 1;
g1 = 3;
g2 = 1;
x=[-2:.01:2];
G1 = gaussian(x,s1,g1);
G2 = gaussian(x,s2,g2);

figure(1)
    plot(x,G1,'--',x,-G2,'--',x,G1-G2)
    xlabel('x')
    box off
    legend('boxoff')
    legend('g_1(x)','g_2(x)','f(x)=g_1(x)+g_2(x)')

end % of main

%%
function out = gaussian(x,s,g)
    out = g*s*exp(-s^2*x.^2)/2;
end
