
function neutropenia_2
% this uses the Matlab routine dde23 to solve delay differential equations
%  for the extended cyclic neutropenia model
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global A beta X 

% parameters
d=2;
X = 10;
beta = 0.5;
A = 1;
 
    
formatSpecF = '%6.2f\n';
 
%specify the output points
tspan = [-10:.01:20];

%specify initial data
lags = [d,d+X];; 
 
sol= dde23(@ddeRHS,lags, @history, tspan) 
 
figure(1)
plot(sol.x,sol.y(1,:),'linewidth',2)
xlabel('t','fontsize',16)
 ylabel('N')
 %title(strcat('\gamma = ',sprintf(formatSpecF,gam)),'fontsize',18)
 axis([0 20 0 2])  
 
set(gca,'linewidth',2.0)
box off
 
 
end

 %%%%%%%%%%%%%%%%%%
 function out = ddeRHS(t,Y,Z)
global  beta X
N=Y(1);

Nlagd = Z(1,1);
NlagdX = Z(1,2);
 
dNdt = -F(NlagdX)*exp(-beta*X) +F(Nlagd)-beta *N;
 
out=[ dNdt];

 end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fout = F(x)
global A 

Fout = A/(1+x^7);
end

function s = history(t)
  s =  [0.8];
end
 