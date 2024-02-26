
%  -------------------------------------------------------------------
%
%   Reponse of the RyR model of Stern to a step increase in calcium.
%
%   For Chapter 7, Section 7.4.4 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 


%% symbolic calculation and plotting of the steady state

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 

syms k1 km1 k2 km2 R O I RI c
global k1val km1val k2val km2val

RI = 1 - R - O - I;
eq1 = km1*O + km2*RI - R*(k1*c*c + k2*c);
eq2 = km2*I + k1*c*c*R - O*(k2*c + km1);
eq3 = k2*c*O + k1*c*c*RI - I*(km2 + km1);

k1val = 35;
km1val = 0.06;
k2val = 0.5;
km2val = 0.01;

[R,O,I] = solve([eq1,eq2,eq3],[R,O,I]);
Ost(c) = subs(O,[k1,km1,k2,km2],[k1val,km1val,k2val,km2val]);
Rst(c) = subs(R,[k1,km1,k2,km2],[k1val,km1val,k2val,km2val]);
Ist(c) = subs(I,[k1,km1,k2,km2],[k1val,km1val,k2val,km2val]);

figure(2)
fplot(Rst,[0.01 0.5],'r','LineWidth',2)
hold on
fplot(Ost,[0.01 0.5],'b','LineWidth',2)
fplot(Ist,[0.01 0.5],'g','LineWidth',2)

% for external plotting
% keep = fplot(Ost,[0.01 0.5]);
% writematrix([keep.XData' keep.YData'],'test1.dat')

%% numerical solutions

%clear all

init = [1,0,0];                 % Assume that c = 0 initially. Probably not very accurate.
tspan = linspace(0,1.5,1000);

figure(1)
c = 1;
[T,Y] = ode15s(@(t,x)rhs(t,x,c),tspan,init);
plot(T,Y(:,2))
%writematrix( [T Y(:,2)],'test3.dat')   % for external plotting
hold on

c = 2;
[T,Y] = ode15s(@(t,x)rhs(t,x,c),tspan,init);
%writematrix( [T Y(:,2)],'test4.dat')   % for external plotting
plot(T,Y(:,2))

c = 3;
[T,Y] = ode15s(@(t,x)rhs(t,x,c),tspan,init);
%writematrix( [T Y(:,2)],'test5.dat')   % for external plotting
plot(T,Y(:,2))
hold off

% calculate the steady-state and peak responses

n = 50;
cc = linspace(0.01,0.5,n);
tspan = linspace(0,50,100000);
for i = 1:n
    c = cc(i);
    [T,Y] = ode15s(@(t,x)rhs(t,x,c),tspan,init);
    peak(i) = max(Y(:,2));
end

figure(2)
plot(cc,peak,'k')
 
% for external plotting
% writematrix([cc',peak'],'test2.dat')



%%
function out = rhs(t,x,c)
global k1val km1val k2val km2val
R = x(1);
O = x(2);
I = x(3);
RI = 1 - R - O - I;

k1 = k1val;
km1 = km1val;
k2 = k2val;
km2 = km2val;

out(1) = km1*O + km2*RI - R*(k1*c*c + k2*c);
out(2) = km2*I + k1*c*c*R - O*(k2*c + km1);
out(3) = k2*c*O + k1*c*c*RI - I*(km2 + km1);
out = out';

end
