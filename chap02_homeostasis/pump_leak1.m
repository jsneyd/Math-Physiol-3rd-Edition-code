%this file makes plots for the original pump-leak model
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
%set parameter values
Ne = 437;
Ke = 20;
gamma = 0.11;
zx = -1;
RTbyF = 25.8;

p = [0.01:.01:14];  %pump rate values could go to 14.3

alpha = (Ne*exp(-3*p)+Ke*exp(2*p*gamma))/(Ne+Ke);
 
 
mu=(1 + sqrt(-alpha*zx^2 + zx^2 + alpha))./(2*(1 - alpha))
y = (-zx+sqrt(zx^2+4*alpha.*mu.^2))./(2*alpha.*mu);


v= -RTbyF*log(y);
 vna = RTbyF*(3*p-log(y));
 vk = RTbyF*(-2*p*gamma-log(y));
figure(1)
plot(p,mu,'r','linewidth',2)
axis([0 14 0 5])
xlabel('Pump rate','fontsize',20)
ylabel('Cell Volume','fontsize',20)
 
figure(2)
plot(p,v,'r',p,vna,'b',p,vk,'g','linewidth',2)
legend('V','V_{Na}','V_K','fontsize',18)
xlabel('Pump rate','fontsize',20)
ylabel('Potential (mV)','fontsize',20)
axis([0 14 -100 100])
 