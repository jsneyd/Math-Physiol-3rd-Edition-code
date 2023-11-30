
function neutropenia_2
% this uses the Matlab routine dde23 to solve delay differential equations
%  for extended cyclic neutropenia model
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global n tauR  b0 K1 gam p0 K2 m tauN A alpha

% parameters
gam = 0.04;
n = 1;
tauR= 2.8;
 K1 = 0.095;
b0 = 8;
p0 = 0.8;
K2 = 0.36;
m = 2;
tauN = 3.5;
A = 10.2;
alpha = 2.4;
 
    
formatSpecF = '%6.2f\n';
 
%specify the output points
tspan = [0:.005:50];

%specify initial data
lags = [tauR,tauN];; 
 
sol= dde23(@ddeRHS,lags, @history, tspan) 
 size(sol)
figure(1)
plot(sol.x,sol.y(1,:),'r' ,sol.x,sol.y(2,:),'b' ,'linewidth',2)
xlabel('t','fontsize',16)
legend('boxoff')
legend('R','N')
 title(strcat('\gamma = ',sprintf(formatSpecF,gam),'fontsize',18)
   
 
set(gca,'linewidth',2.0)
box off
 
figure(2)
 plot(sol.y(1,:), sol.y(2,:))
 xlabel('R')
 ylabel('N')

end

 %%%%%%%%%%%%%%%%%%
 function out = ddeRHS(t,Y,Z)
global  gam A alpha tauR
R=Y(1);
N = Y(2);
RlagR = Z(1,1);
NlagN = Z(2,2);
RlagN = Z(1,2);

dRdt = 2*exp(-gam*tauR)*beta(RlagR)*RlagR - (phi(N)+beta(R))*R;
dNdt = A*phi(NlagN)*RlagN-alpha *N;
 
out=[dRdt;dNdt];

 end 
function gbout = beta(x)
global K1 n  b0

gbout = b0* K1^n/(K1^n+x^n);
end

function gpout = phi(x)
global K2 m  p0

gpout = p0* K2^m/(K2^m+x^m);
end

function s = history(t)
  s =  [2.8;3];2*ones(2,1);
end
 