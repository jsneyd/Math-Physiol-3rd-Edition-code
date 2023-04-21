clear all
close all
clc

global p

getparams
init = [p.h0,0,8];

% Staircase experiment
G = [0,50,100,150,200];     % glucose levels  
times = [0 2 7 12 17];     % times of step application
tspan = linspace(0,22,500);
[T,Y] = ode15s(@(t,x)rhs(t,x,G,times),tspan,init);
F = Y(:,p.ng+1);
I = Y(:,p.ng+2);
figure(1)
plot(T,F)
%save('pedersen1.mat','T','F')

%same stimulation twice
G = [0,300,0,300];     % glucose levels  
times = [0 5 60 65];     % times of step application
tspan = linspace(0,75,500);
[T,Y] = ode15s(@(t,x)rhs(t,x,G,times),tspan,init);
F = Y(:,p.ng+1);
I = Y(:,p.ng+2);
figure(2)
plot(T,F)
%save('pedersen2.mat','T','F')






%%
function out = rhs(t,x,G,times)
global p
h = x(1:p.ng)';
F = x(p.ng+1);
I = x(p.ng+2);

dum = find(times <= t);
p.glucose = G(dum(end));

M = p.M1*p.glucose^p.m/(p.KM^p.m + p.glucose^p.m) + p.M0;

int = trapz(p.g(p.g<p.glucose),h(p.g<p.glucose));
intinf = trapz(p.g,h);
out(1:p.ng) = p.pplus*I*p.phi - p.pminus*h - p.fplus*h.*heaviside(p.glucose-p.g);
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