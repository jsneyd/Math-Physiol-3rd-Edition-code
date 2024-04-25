%-------------------------------------------------------------------

% Matlab code for simulating the Goldbeter-Lefever model of glycolytic oscillations.

% For Chapter 1 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------


function GL_glycolysis_sim
close all
clear all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
 

% this code is to simulate the Goldbeter-Lefever enzyme model

%parameters
 
 par.nu=190;
 par.eta=120;
 
 
% set up the integration
tspan = [0:0.0005:1];            % time interval for the solution
IC = [37,1];                        % initial condition
[t,U] = ode15s(@(t,y)rhs(t,y,par),tspan,IC);

figure(1)
plot(t,U(:,1),'r',t,U(:,2),'b--','linewidth',2)
legend('boxoff')
legend('\sigma_1','\sigma_2','fontsize',18,'location','northwest')
xlabel('time','fontsize',18)
ylabel('concentration','fontsize',18)

figure(2)

% null clines
w=[0.02:0.1:20];


u1=par.nu./(1+w).^2;

u2 =par.eta*w./(1+w).^2;
 
plot(U(:,1),U(:,2),'r',u2,w,'b--',u1,w,'g--','linewidth',2)
xlabel('\sigma_1','fontsize',20)
ylabel('\sigma_2','fontsize',20)
legend('boxoff')
legend('solution','\sigma_1 nullcline','\sigma_2 nullcline','fontsize',18)
axis([0 40 0 20])
hold off

%HB curve
eta = ones(201,1)*[0: 1:200];
nu = [0:1:200]'*ones(1,201);

F=-eta.^4 + eta.^3.*nu - eta.^3 - 3*eta.^2.*nu - 3*eta.*nu.^2 - nu.^3;
figure(3)
contour(nu,eta,F,[0 0],'linewidth',2)
hold on
plot(par.nu,par.eta,'*')
xlabel('\nu')
ylabel('\eta')
text(50,150,'stable','fontsize',20)
text(100,60,'unstable','fontsize',20)
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=rhs(t,y,par)
 
u=y(1);
w=y(2);
f =u*(1+w)^2;
Fu = par.nu-f;
Fw =  f -par.eta*w;
 
out = [Fu Fw]';
