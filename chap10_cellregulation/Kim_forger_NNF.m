
% Code to solve the NNF circadian clock model of Kim and Fprger.
% Used for Keener and Sneyd, Mathematical Physiology, Third Edition,
% Chapter 10, Exercise ??


clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

global Kd alpha1 alpha2 alpha3 bet1 bet2 bet3  AT del
%specify parameter values
 formatSpecF = '%6.2f\n';
alpha1 = 20;
alpha2 = 1;
alpha3 = 1;
bet1 = 1
bet2 = 0;
bet3 = 1;
 
 Kd=1.e-2;
 AT= 1;
   
 del = 0.5
  
%specify the output points
tspan = [-100:.1:30];

%specify initial data

y0 = [0.5;0.5;0.5;1;0.5];
 
[T,Y] = ode23(@deRHS,tspan, y0);  
figure(1)
plot(T,Y(:,1), T,Y(:,2), T, Y(:,3),'linewidth',2)
 title(strcat('K_1 = ',sprintf(formatSpecF,Kd),', \alpha_1 = ',sprintf(formatSpecF,alpha1)),'fontsize',18)

xlabel('time','fontsize',16)
legend('boxoff')
legend('M','P_c','P_T')
%text(16,3,'PER','fontsize',16)
%text(11,1.5,'mRNA','fontsize',16)
axis([0 20 5 25])
set(gca,'linewidth',1.5)
box off
%legend('mRNA','Per')

 

%%
function F_prime=deRHS(t,y)
global Kd alpha1 alpha2 alpha3 bet1 bet2 bet3  AT del

M = y(1);
Pc = y(2);
PT = y(3);
AT = y(4);
V = y(5);
 

A = (AT-PT-Kd+((AT-PT-Kd)^2+4*AT*Kd)^0.5)/2;
Mp=alpha1*A-bet1*M;
Pcp=alpha2*M-(bet2+alpha3)*Pc;
PTp=alpha3*Pc-bet3*PT;
 
ATp = del*(1/(1+V) - A)
Vp = del*(A - V)

 
F_prime = [Mp;Pcp;PTp;ATp;Vp];
end
