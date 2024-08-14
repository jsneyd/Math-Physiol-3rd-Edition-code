
%    -------------------------------------------------------------------
%
%     Model of unregulated glomerular filtration.
%
%     For Chapter 17, Section 17.1 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function glomerulus2

clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

global P1 Ra alpha KfL pii Pa

% use typical parameters to find Ra, Re, Rd and KfL
P1 = 60;
P2 = 18;
Pa = 100;
Pe = 0;
Pd = 18;
pii = 25; %mm Hg
Qi = 650;
Qd = 125;
Qe = Qi-Qd;
alpha = pii/(P1-P2);

Ra = (Pa-P1)/Qi;
Re = (P1-Pe)/Qe;
Rd = (P2-Pd)/Qd;
KfL = -(Qe/Qi + alpha*log((Qe/Qi-alpha)/(1-alpha)) -1) *(alpha *Qi)/pii;

% now allow Pa to vary and find the corresponding values of Qe and Qd
n = 50;
Pa_values = linspace(67,160,n);
options = optimset('Display','off');
for i=1:n
    Pa = Pa_values(i);
    Qe(i) = fsolve(@(Qe)rhs_root(Qe),100,options);
end

Qi = (Pa_values - P1)/Ra;
plot(Pa_values,0.1*Qe,'r',Pa_values,Qi-Qe,'b','LineWidth',2)
xlabel('P_a')
legend('Q_e/10','Q_d')

end % of main


%%
function out = rhs_root(Qe)
    global P1 Ra alpha KfL pii Pa
    Qi = (Pa - P1)/Ra;
    out = Qe/Qi + alpha*log((Qe/Qi-alpha)/(1-alpha)) - 1 + KfL*pii/(alpha*Qi);
end
