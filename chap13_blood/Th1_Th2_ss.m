% this plots the Th1 Th2 steady state curves for Fig. 13.17
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);


% set parameters from Table 13.5

sig = 5;
kap = 1;
rho = 1;
beta = 0.05;
alpha = 5;
g1 = 1;
g2 = 0.5;
mu = 3;
g = 2;
n = 4;

x1=[0.287:.001:5];
% select a value of S1
S1 = 0;

f1 = (-beta + mu*x1)./(alpha*x1.^n./(kap^n+x1.^n) +sig*S1/(rho+S1)) ;

x2=(g2./f1)-g2;

f2=(-beta + mu*x2).*(g1+x1)/g1-alpha.*x2.^n./(kap^n+x2.^n);

Sig2=f2;
S2=rho* Sig2./(sig-Sig2);

figure(1)
plot(S2,x1,'--',S2,x2)
legend('boxoff')
legend('T-bet','GATA-3')
axis([0 2 0 2.5])
box off
xlabel('S_2')

%figure(2)
%plot(x1,x2)

