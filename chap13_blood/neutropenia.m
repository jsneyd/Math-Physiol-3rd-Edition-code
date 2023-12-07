
function neutropenia
% this uses the Matlab routine dde23 to solve delay differential equations
%  for cyclic neutropenia
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global n tau del b0 K1 gam

% parameters
gam = 0.33;
n = 4;
tau= 1.39;
 del = 0.08;
 K1 = 3.07;
b0 = 1.62;
 
 
%specify the output points
tspan = [0:.01:200];

%specify initial data
lags = tau; 
 
sol= dde23(@ddeRHS,lags, @history, tspan) 
 size(sol)
figure(1)
plot(sol.x,sol.y(1,:),'r' ,'linewidth',2)
xlabel('t','fontsize',16)
ylabel('R')
 
 
set(gca,'linewidth',1.5)
box off
 
 
end

 %%%%%%%%%%%%%%%%%%
 function out = ddeRHS(t,Y,Z)
 global   n tau del b0 K1 gam 
R=Y(1);

Rlag = Z(1);

dRdt = 2*gs(Rlag)*Rlag*exp(-gam*tau) - (del+gs(R))*R;

 
out=[dRdt];

 end 
function gsout = gs(x)
global K1 n  b0

gsout = b0* K1^n/(K1^n+x^n);
end

function s = history(t)
  s = 4.7;%ones(2,1);
end
 