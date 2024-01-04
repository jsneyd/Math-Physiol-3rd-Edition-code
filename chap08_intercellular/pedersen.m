%  -------------------------------------------------------------------
%
%   Code for the Pedersen model of secretory granule exocytosis.
%
%   For Chapter 8, Section 8.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

global p

getparams
init = [p.h0,0,8];  % the first ng components are h. Then F, then I.

% Staircase experiment
p.stim_choose = 1;  % 1 for the staircase, 2 for the repeated stimulation
tspan = linspace(0,21,1000);
[T,Y] = ode15s(@(t,x)rhs(t,x),tspan,init);
F = Y(:,p.ng+1);
I = Y(:,p.ng+2);
h = Y(:,1:p.ng);
intinf = trapz(p.g,h');

figure(1)
plot(T,F)
ylabel('insulin secretion (\mug/min)')
yyaxis ('right')
plot(T,glucose(T),'--')
ylabel('glucose (mg/100 ml)')
xlabel('time(min)')
xlim([0,21])
legend('boxoff')
legend('model','glucose','Location','Northwest')

%figure(6)
%plot(T,intinf)

%save('pedersen1.mat','T','F')

%same stimulation twice
p.stim_choose = 2;  % 1 for the staircase, 2 for the repeated stimulation
tspan = linspace(0,75,500);
[T,Y] = ode15s(@(t,x)rhs(t,x),tspan,init);
F = Y(:,p.ng+1);
I = Y(:,p.ng+2);

figure(2)
plot(T,F)
ylabel('insulin secretion (\mug/min)')
yyaxis ('right')
plot(T,glucose(T),'--')
ylabel('glucose (mg/100 ml)')
xlabel('time(min)')

%save('pedersen2.mat','T','F')

%% glucose stimulation as a function of time
function out=glucose(t)
global p
if p.stim_choose==1
    G = [0,50,100 150 200];     % glucose levels  
    times = [0 2 7 12 17];     % times of step application 
    out = G(1)*(heaviside(t-times(1)) - heaviside(t-times(2))) + ...
          G(2)*(heaviside(t-times(2)) - heaviside(t-times(3))) + ...
          G(3)*(heaviside(t-times(3)) - heaviside(t-times(4))) + ...
          G(4)*(heaviside(t-times(4)) - heaviside(t-times(5))) + ...
          G(5)*(heaviside(t-times(5)));
end

if p.stim_choose==2
    G = [0,300,0,300];     % glucose levels  
    times = [0 5 60 65];     % times of step application
    out = G(1)*(heaviside(t-times(1)) - heaviside(t-times(2))) + ...
      G(2)*(heaviside(t-times(2)) - heaviside(t-times(3))) + ...
      G(3)*(heaviside(t-times(3)) - heaviside(t-times(4))) + ...
      G(4)*heaviside(t-times(4));
end

end


%%
function out = rhs(t,x)
global p
h = x(1:p.ng)';
F = x(p.ng+1);
I = x(p.ng+2);

p.glucose = glucose(t);
M = p.M1*p.glucose^p.m/(p.KM^p.m + p.glucose^p.m) + p.M0;

int = trapz(p.g(p.g<p.glucose),h(p.g<p.glucose));
intinf = trapz(p.g,h);
out(1:p.ng) = p.pplus*I*p.phi - p.pminus*h - p.fplus*h.*(p.glucose>p.g);
out(p.ng+1) = p.fplus*int - p.k*F;
out(p.ng+2) = M - p.r*I - p.pplus*I + p.pminus*intinf;

out = out';

end

%%
function getparams
global p

p.Xmax = 1.65;
p.n = 3.3;
p.K = 150;

p.M1 = 1.25;
p.m = 12;
p.KM = 185;
p.M0 = 0.02;

p.fplus = 6.2;
p.k = 1.09;
p.pplus = 0.021;
p.pminus = 0.11;
p.r = 0.0023;

p.ng = 600;  % Number of values of g
p.g = linspace(0,600,p.ng);
p.phi = p.n*p.g.^(p.n-1)*p.K^p.n./((p.K^p.n + p.g.^p.n).^2);
p.h0 = p.Xmax*p.phi;
%p.h0 = linspace(0,0,p.ng);   % This also works
end

