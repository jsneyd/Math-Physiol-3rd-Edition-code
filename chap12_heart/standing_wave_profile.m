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
mu = exp(-L/lambda_g)
m=mu;
ratc=qi/Q; 
E = exp(sQ);
LAM = 2*(mu-E)*(mu-1/E)/(mu*(E-1/E));

A=-E/(E^2*LAM + 2*E^2 - 4*E*m - LAM + 2)
 
B=A;
C= -(2*E^2*LAM*m^2 + 2*E^2*m^2*ratc - 2*E^2*LAM*m + 2*E^2*m^2 - 4*E*m^3 ...
    - 4*E*m^2*ratc - 4*E^2*m - 2*E^2*ratc + 8*E*m^2 - 2*LAM*m^2 ...
    + 2*m^2*ratc + 2*E^2 - 4*E*m + 4*E*ratc + 2*LAM*m + 2*m^2 - 4*m - 2*ratc + 2) ...
    /((E^2*LAM + 2*E^2 - 4*E*m - LAM + 2)*(m^2 - 2*m + 1))
 


x = [0:.05:1]*L;
% for positive n
phi = (mu-1/E)*exp(lambda*x) + (mu - E)*exp(-lambda*x);

Pk = 2*qe/Q*(mu-1/E)*(mu-E)/(mu-1);
phi_i = phi*qi/Q+Pk;
phi_e = -phi*qe/Q + Pk;

% for negative n

phin = (1/mu-1/E)*exp(lambda*x) + (1/mu - E)*exp(-lambda*x);

Pkn = 2*qe/Q*(1/mu-1/E)*(1/mu-E)/(1/mu-1);
phi_in =  (phin*qi/Q+Pkn);
phi_en = -phin*qe/Q + Pkn;

Vi = [];
Ve = [];
y = [];

 
for j =  -25:-1
Vi = [Vi,B*mu^(-j)*phi_in+1+C];
 Ve=[Ve,B*mu^(-j)*phi_en+C];
y=[y,j*L+x];
end



for j = 0:25
    Vi = [Vi,A*mu^j*phi_i];
    Ve = [Ve,A*mu^j*phi_e];
    y = [y,j*L+x];
end

figure(1)
plot(y,Ve,y,Vi,y,Vi-Ve,'linewidth',2);

xlabel('Length (cm)', 'fontsize',18)
ylabel('Potential', 'fontsize',18)
text(-.1,-0.05,'V_e','fontsize',18)
text(-.1,0.6,'V_i','fontsize',18)
text(-0.05,0.8,'V','fontsize',18)
 axis([-0.3 0.3 -0.2 1])
set(gca,'linewidth',2.0)
box off
test=[y' Vi' Ve'];
save('test.dat','-ascii','test')

 
