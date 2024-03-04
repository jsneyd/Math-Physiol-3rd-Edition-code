
%  -------------------------------------------------------------------
%
%   Make the plot for Fig. 13.10.
%
%   For Chapter 13, Section 13.1.4 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
% parameter sets
taulist = [3.8, 1.2];
alist = [0.36, 0.53];
klist = [3.15,4.38];

omegabya = [0:.01:pi/4];


for j = 1:2
    tau = taulist(j);
    a = alist(j);
    kp1 = klist(j);
    omega = a*omegabya;
    arg = omega.*tau+kp1*atan(omegabya);
    alpha = -omega./tan(arg);
    mu = -omega./(cos(atan(omegabya)).^kp1.*sin(arg));
    ndx = find(alpha>0);
    
    figure(1)
    if (j==1)
        plot(alpha(ndx),mu(ndx))
    end
    if (j==2)
        plot(alpha(ndx),mu(ndx),'--')
    end
    hold on

end
axis([0 4 -9 0])
box off

text(0.5,-5,'Unstable','fontsize',18)
text(3,-2,'Stable','fontsize',18)
text(2.7,-5,'Normal Human','fontsize',18)
text(2.55,-8,'CN Human','fontsize',18)

plot([1.75,1.75],[-9,0],'k','linewidth',0.5)
plot([2.5,2.5],[-9,0],'k','linewidth',0.5)
xlabel('\alpha (day^{-1})')
ylabel('\mu (day^{-1})')
hold off
