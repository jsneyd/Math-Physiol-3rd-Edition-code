
%  -------------------------------------------------------------------
%
%   Code to simulate the SNF circadian clock model of Kim and Forger.
%
%   For Chapter 10, Section 10.3.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function Kim_forger_SNF

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

global Kd alpha1 alpha2 alpha3 bet1 bet2 bet3  AT
%specify parameter values
 formatSpecF = '%6.3f\n';
alpha1 = 20;
alpha2 = 1;
alpha3 = 1;
bet1 = 1
bet2 = 0;
bet3 = 1;

 Kd=1.e-3;
 AT= 1;

%specify the output points
tspan = [-24:.1:30];

%specify initial data

y0 = [0.5;0.5;0.5];

[T,Y] = ode23(@deRHS,tspan, y0);
figure(1)
plot(T,Y(:,1), T,Y(:,2),T,Y(:,3),'linewidth',2)
 title(strcat('K_1 = ',sprintf(formatSpecF,Kd),', \alpha_1 = ',sprintf(formatSpecF,alpha1)),'fontsize',18)


%xlabel('time (hours)','fontsize',16)
%text(16,3,'PER','fontsize',16)
%text(11,1.5,'mRNA','fontsize',16)
 axis([0 20 0 3])
 legend('boxoff')
 legend('M','P_c','P_T')
 xlabel('t','fontsize',18)
set(gca,'linewidth',2)
box off



% figure(2)
% loglog(P,M1,'b--',P,M2,'r--',Y(:,2),Y(:,1),'linewidth',2)
% ylabel('mRNA','fontsize',16)
% xlabel('Per','fontsize',16)
% text(2.5,.1,'dM/dt = 0','fontsize',16)
% text(2,4,'dP/dt = 0','fontsize',16)
% %axis([.001 100 .01 100])
% set(gca,'linewidth',1.5)
% box off

end % of main


%%
function F_prime=deRHS(t,y)
global Kd alpha1 alpha2 alpha3 bet1 bet2 bet3  AT

M = y(1);
Pc = y(2);
P = y(3);


A = (AT-P-Kd+((AT-P-Kd)^2+4*AT*Kd)^0.5)/2;
Mp=alpha1*A-bet1*M;
Pcp=alpha2*M-(bet2+alpha3)*Pc;
Pp=alpha3*Pc-bet3*P;


F_prime = [Mp;Pcp;Pp];
end
