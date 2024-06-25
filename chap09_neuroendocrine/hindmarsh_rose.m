
%  -------------------------------------------------------------------
%
%   Code for the Hindmarsh-Rose model of electrical bursting.
%
%   For Chapter 9, Section 9.1.3 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  -------------------------------------------------------------------

close all
clear all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);

global   Iapp  r s  x1

% parameters
r=0.001;
s=4;
x1 = -(1+sqrt(5))/2;

zstart=[0.2,1.8,2];         % initial data
total=1000;
tstep = 0.01;
Ilist=[0.4,2,4];
tspan = [0:tstep:total];    % specify the output points

for j=1:length(Ilist)
    Iapp=Ilist(j)  
    init=[0,0,zstart(j)];
    [T,S] = ode15s(@deRHS,tspan, init, odeset('maxstep',10));  
    
    figure(2*j-1)
    plot(T ,S(:,1),T,S(:,3) ,'--')
    formatSpecF = '%6.2f\n';
     
    title(strcat('I_{app} = ',sprintf(formatSpecF,Iapp)),'fontsize',18)
    legend('boxoff')
    legend('x','z')
    ylabel('x')
    xlabel('T')
    box off

    figure(2*j)
    plot(S(:,3),S(:,1) )
    box off
    formatSpecF = '%6.2f\n';
     
    title(strcat('I_{app} = ',sprintf(formatSpecF,Iapp)),'fontsize',18)
    xlabel('z')
    ylabel('x')
end

%%
function s_prime=deRHS(t,sol) 
global  Iapp  r s x1
x=sol(1); 
y=sol(2);
z=sol(3);

zp = r*(s*(x-x1)-z);
xp = y - x^3 + 3*x^2 + Iapp -z  ;
yp = 1 - 5*x^2 - y;
s_prime =[xp;yp;zp];
end