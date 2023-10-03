%this  file calculates the potentials in a passive cardiac cable.

%The sawtooth potential
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)

L = 0.012;
qi = 5.47e-3;
qe = 0.5 *qi;
Q = qi+qe;
sQ = sqrt(Q);
lambda = sQ/L;
lambda_g = 0.09;
mu = exp(-L/lambda_g);



E = exp(sQ);

A=mu/(2*(1-mu^2));
B=mu^2*A;
C=-2*qe/(qe+qi)*A*(E-1)^2/E;
% for negative n

x = [0:.05:1]*L;
phin = (1/mu-1/E)*exp(lambda*x) + (1/mu - E)*exp(-lambda*x);

Pkn = 2*qe*(1/mu-1/E)*(1/mu-E)/(1/mu-1);
phi_in = phin*qi+Pkn;
phi_en = -phin*qe + Pkn;

Vi = [];
Ve = [];
y = [];

figure(4)
plot(x,phin,x,phi_in,x,phi_en)

for j =  -15:-1
Vi = [Vi,B*mu^(-j)*phi_in+1+C];
 Ve=[Ve,B*mu^(-j)*phi_en+C];
y=[y,j*L+x];
end

phi = (mu-1/E)*exp(lambda*x) + (mu - E)*exp(-lambda*x);

Pk = 2*qe*(mu-1/E)*(mu-E)/(mu-1);
phi_i = phi*qi/(qi+qe)+Pk;
phi_e = -phi*qe/(qi+qe) + Pk;


for j = 0:15
    Vi = [Vi,A*mu^j*phi_i];
    Ve = [Ve,A*mu^j*phi_e];
    y = [y,j*L+x];
end

figure(1)
plot(y,Ve,y,Vi,y,Vi-Ve,'linewidth',2);

xlabel('Length (cm)', 'fontsize',18)
ylabel('Potential', 'fontsize',18)
text(.05,0,'V_e','fontsize',18)
text(.05,-.8,'V_i','fontsize',18)
text(.05,-1.4,'V','fontsize',18)
 
set(gca,'linewidth',2.0)
box off
test=[y' Vi' Ve'];
save('test.dat','-ascii','test')



