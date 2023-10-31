clear all
close all
clc

L = 3;   % Units of cm
l = 0.0035;
N = 1000;
x = linspace(0,L,N);
w = 800;
sigma = l/L;

k0 = 1e7;
r0 = 3000;
F0 = 1000;
lam = 1.5;
rho = 1;
k = k0*exp(-lam*x); 
r = r0*exp(-lam*x);
Z = (r + k/(1i*w));
Z0 = (r0 + k0/(1i*w));
Y = 2./Z;
eta0 = F0/Z0;

intY = 2*1i*w*(exp(lam*x) - 1)/(lam*(1i*w*r0 + k0));
eta = -1i*rho*eta0.*Y.*exp(-w*rho*intY);


figure(1)
plot(x,eta,'r','LineWidth',1)
hold on
plot(x,abs(eta),'b','LineWidth',2)
hold off

% figure(2)
% % test with the analytic solution
% beta = 2*rho*(w^2)*(w*r0 + 1i*k0)/(lam*(k0^2 + w^2*r0^2));
% xi = w*rho/(1i*w*r0+k0);
% eta2 = 2*eta0*xi*exp(lam*x + beta*(1-exp(lam*x)));
% plot(x,eta2,'r','LineWidth',2)
% hold on
% plot(x,abs(eta2),'b','LineWidth',2)



