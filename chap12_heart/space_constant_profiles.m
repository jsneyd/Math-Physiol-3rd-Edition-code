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

x = [0:.05:1]*L;

phi = (mu-1/E)*exp(lambda*x) + (mu - E)*exp(-lambda*x);

Pk = 2*qe/Q*(mu-1/E)*(mu-E)/(mu-1);
phi_i = phi*qi/Q+Pk;
phi_e = -phi*qe/Q + Pk;
Vi = [];
Ve = [];
y = [];

sc =1.e2;
for j = 0:15;
    Vi = [Vi,sc*mu^j*phi_i];
    Ve = [Ve,sc*mu^j*phi_e];
    y = [y,j*L+x];
end

figure(1)
plot(y,Ve,y,Vi,y,Vi-Ve,'linewidth',2);
 
xlabel('Length (cm)', 'fontsize',18)
ylabel('Potential', 'fontsize',18)
text(.05,0,'V_e','fontsize',18)
text(.05,-18,'V_i','fontsize',18)
text(.05,-11,'V','fontsize',18)
 
set(gca,'linewidth',2.0)
box off
test=[y' Vi' Ve'];
save('test.dat','-ascii','test')



