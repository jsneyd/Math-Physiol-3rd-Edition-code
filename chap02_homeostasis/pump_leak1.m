%-------------------------------------------------------------------

% Matlab code for plotting the solutions of the pump-leak model.

% For Chapter 2, Fig. 2.17 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------
clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 
  
%set parameter values
Ne = 437;
Ke = 20;
gamma = 0.11;
zx = -1;
RTbyF = 25.8;

p = [0.01:.01:14];  %pump rate values could go to 14.3

alpha = (Ne*exp(-3*p)+Ke*exp(2*p*gamma))/(Ne+Ke);
 
 
mu=(1 + sqrt(-alpha*zx^2 + zx^2 + alpha))./(2*(1 - alpha));
y = (-zx+sqrt(zx^2+4*alpha.*mu.^2))./(2*alpha.*mu);


v= -RTbyF*log(y);
 vna = RTbyF*(3*p-log(y));
 vk = RTbyF*(-2*p*gamma-log(y));
figure(1)
plot(p,mu,'r','linewidth',2)
axis([0 14 0 5])
xlabel('Pump rate, P','fontsize',20)
ylabel('Cell Volume, \mu','fontsize',20)
box off

figure(2)
plot(p,v,'r',p,vna,'b',p,vk,'g','linewidth',2)
legend('boxoff')
legend('V','V_{Na}','V_K','fontsize',18)
xlabel('Pump rate, P','fontsize',20)
ylabel('Potential (mV)','fontsize',20)
axis([0 14 -100 100])
 box off



% Now the modified model with P = rho u^3

u = y.*exp(-3*p);

rho=p./ u.^3;
% plot everything as a function of rho
 
figure(3)
plot(rho,mu,'r','linewidth',2)
axis([0 4  0 5])
xlabel('Pump rate, \rho','fontsize',20)
ylabel('Cell Volume, \mu','fontsize',20)
box off
 

 figure(4)
plot(rho,v,'r',rho,vna,'b',rho,vk,'g','linewidth',2)
legend('boxoff')
legend('V','V_{Na}','V_K','fontsize',18,'location', 'best')
xlabel('Pump rate, \rho','fontsize',20)
ylabel('Potential (mV)','fontsize',20)
 axis([0 400 -80 60])
 box off
  

 exercise_on = 0;
% Now the modified model with P = rho u^3/y for exercise 2.30
if (exercise_on==1)
 
rho=y.*p./ u.^3;
% plot everything as a function of rho
 
figure(5)
plot(rho,mu,'r','linewidth',2)
axis([0 4  0 5])
xlabel('Pump rate','fontsize',20)
ylabel('Cell Volume','fontsize',20)
box off

 figure(6)
plot(rho,v,'r',rho,vna,'b',rho,vk,'g','linewidth',2)
legend('boxoff')
legend('V','V_{Na}','V_K','fontsize',18,'location', 'best')
xlabel('Pump rate \rho','fontsize',20)
ylabel('Potential (mV)','fontsize',20)
 axis([0 400 -80 60])
 box off

end

 
