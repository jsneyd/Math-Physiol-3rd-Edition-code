%  -------------------------------------------------------------------
%
%   The feedforward motif.
%
%   For Chapter 10, Section 10.1.3 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function motif_FF

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
init=[0 0]
tspan=linspace(0,50,5000);

[T,Y] = ode15s(@rhs,tspan,init);
R = Y(:,1);
X = Y(:,2);
S = (T>10) + (T>20) +  (T>30);

plot(T,R)
ylabel('R')
yyaxis('right')
plot(T,S)
ylabel('S')
xlabel('t')
legend('boxoff')
legend('R','S','location','southeast')
set(gca,'linewidth',1.5)
box off
out = [T R S];
%save('test.dat','out','-ascii')

end % of main


%%
function out = rhs(t,Y)
R = Y(1);
X = Y(2);

n=2;
K=2;
S =  (t>10) +  (t>20) + (t>30);
phi = K*X^n/(K+X^n);
%phi = X; % for exercise 1.9
out(1) = S - phi*R;
out(2) = S - X;
out = out';

end
