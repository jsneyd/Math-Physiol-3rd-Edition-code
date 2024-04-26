
%    -------------------------------------------------------------------
%
%     Co-current/countercurrent comparison.
%
%     For Chapter 17, Section 17.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

clear all; close all; clc;
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0); 


alpha = 10;
x = [0:.01:1];

Cco = (1-exp(-2*alpha*x))/2;
Ccount = alpha*x/(1+alpha);

figure(1)
plot(x,Cco,x,Ccount)
xlabel('x/L')
ylabel('C_1/C')
text(.2,.55,'cocurrent','fontsize',18)
text(.55,.8,'countercurrent','fontsize',18)

