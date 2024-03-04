

%  -------------------------------------------------------------------
%
% Code to solve red blood cell delay differential equation
% using method of lines.
%
%   For Chapter 13, Section 13.1.3 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------


function delay_de
global X d A dx dy M N n0 N0

close all
%clear all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 

 
% first make some  red blood cell production plots

% steady state analysis
N = [0:.01:3];
F = 1./(1+N.^7);
blist = [0.2,0.5,0.8];

figure(1)
    plot(N,F,'linewidth',2)
    xlabel('U','fontsize',20)
    ylabel('Production rate f(U)/A','fontsize',20)
    hold on
    for j = 1:3
        b = blist(j);
        plot(N,b*N,'linewidth',2)
    end
    axis([0 2 0 1])
    text(1.6,0.28,'\beta = 0.2','fontsize',20)
    text(1.6,0.75,'\beta = 0.5','fontsize',20)
    text(1.04,0.8,'\beta = 0.8','fontsize',20)
    box off
    hold off
 
figure(2)
    plot(N./F,N,'linewidth',2)
    ylabel('U_0','fontsize',20)
    xlabel('XA','fontsize',20)
    axis([0 4 0 1.4])
 
% There are three parameters
d = 7;  % time delay
X = 50; %lifetime of red blood cells
A = 1/15;  %   value for f(0);
A = 0.1;

% plot 3: stability diagram
x = [0:.01:16];
f = pi*x./((2+x).*sin(2*pi./(2+x)));
p = 7;
N0p = f./(p-f);
N0 = N0p.^(1/p);
dA = N0.*(1+N0p)./x;
np = find(N0p>0);

figure(3)
    plot(x(np),dA(np) ,'linewidth',2)
    hold on
    % plot( X/d,d*A,'*','linewidth',5)
    hold off
    xlabel('X/d','fontsize',20)
    ylabel('d A','fontsize',20)
    axis([0 14 0 2])
    text(5,1.2,'Unstable','fontsize',20)
    text(5,0.2,'Stable','fontsize',20)
    box off
  
% now do the time dependent simulation
N = 100;  % number of grid points for n
M = 10; %Number of grid points for N
dx  = X/(N-1); %grid spacing for n
dy = d/(M-1); % grid spacing for N
 
tstep = 1;
t_end = 500;
tspan = [0:tstep:t_end];
 
s0 = [ones(M,1);ones(N,1)*dx/X];  %initial data   
[T,S] = ode15s(@deRHS,tspan,s0);

n0 = A./(1+S(:,M).^7);
n = [n0,S(:,M+1:M+N)];
N0 = (sum(n(:,2:end-1),2)+n0/2+n(:,N+1)/2)*dx;  %The integral
figure(4)
    plot(T,N0, 'linewidth',2)
    xlabel('time (days)','fontsize',20)
    ylabel('N(t)','fontsize',20)
    axis([0 500 0 2.5]) 

end % of main

%%
function s_prime=deRHS(t,s)
global  X d A dx dy M N n0 N0
 
    n0 = A/(1+s(M)^7);
    n = [n0;s(M+1:M+N)];
    N0 = (sum(n(2:end-1))+n0/2+n(end)/2)*dx;  %The integral
    Nt = [N0;s(1:M)];
    FN = (Nt(1:M)-Nt(2:M+1))/dy;
    Fn =  (n(1:N)-n(2:N+1))/dx;
    
    s_prime = [FN; Fn];
end
 