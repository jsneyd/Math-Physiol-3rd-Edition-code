%this code is to plot the response to a stimulus in a bidomain cable.

%These are the same parameters as used in the space constant calculation
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
sQ = sqrt(Q);  % = lam *L
lam = sQ/L;
lam_g = 0.09;  %the space constant
lam_g = 0.1;

mu = exp(L/lam_g);
Rm = 7000;
Rc = 150;
Ai = 4e-6;
p = 2*sqrt(pi*Ai);
rc = qi*Rm/(p*L^2);
re = qe*Rm/(p*L^2);
J = 1;
I = lam*J/re;

S = sinh(sQ);
C = cosh(sQ);

%calculate Kg
Kg = ((mu^2+1)/mu -2*C)/S;

% pick number of cells
N = 95;

%find c1 and c2
y1 = [S+Kg*C;mu-C];
x1 = [S;mu-C];
y2 = [S+Kg*C;1/mu-C];
x2 = [S;1/mu-C];
Rhs = [1-C;Kg*(1-C)-S];
nm1 = y1'*x1;
nm2 = y2'*x2;
c1 = y1'*Rhs/nm1;
c2 = y2'*Rhs/nm2;


a0 = J*(mu^(-N)-1)*(c1+c2*mu^(1-N))/((mu-1)*(mu^(-2*N)-1));
j = [0:N];
a1 = -a0*mu.^j+J*c1*(mu.^j-1)/(mu-1);
a2 = a0*mu.^(-j)+J*c2*(mu.^(-j)-1)/(1/mu-1);
cfs = x1*a1+x2*a2;
b = cfs(2,:);
a = cfs(1,:)-J;
nx = 10;
x = [0:1/nx:1]*L;
n = (nx+1)*(N+1);

gamma(1) = 0;
for j = 1:N 
    gamma(j+1) = gamma(j) + (sQ*rc*J +re*Kg*(J*(1-C)+C*cfs(1,j)+S*b(j)))/(re+rc);
end


phi = sinh(lam*x)'*a+cosh(lam*x)'*b;
phi_i = rc*phi/(re+rc) + (rc*J*lam*x'*ones(1,N+1))/(re+rc) + ones(nx+1,1)*gamma;
phi_e = -re*phi/(re+rc) + (rc*J*lam*x'*ones(1,N+1))/(re+rc) + ones(nx+1,1)*gamma;
dphi = lam*cosh(lam*x)'*a+lam*sinh(lam*x)'*b;
dphi_i = rc*dphi/(re+rc) + (rc*J*lam)/(re+rc);
dphi_e = -re*dphi/(re+rc) + (rc*J*lam)/(re+rc);

n = (nx+1)*(N+1);
phi = reshape(phi,1,n);
phi_i = reshape(phi_i,1,n);
phi_e = reshape(phi_e,1,n);
dphi = reshape(dphi,1,n);
dphi_i = reshape(dphi_i,1,n);
dphi_e = reshape(dphi_e,1,n);
current = dphi_i/rc+dphi_e/re;

y = [];
for j = 1:N+1
    y = [y,j*L+x];
end
figure(1)
subplot(2,1,1)
plot(y,phi,'linewidth',2) %y,phi_i,y,phi_e,
xlabel('Length (cm)','fontsize',16)
ylabel('Potential','fontsize',16)
subplot(2,1,2)
plot(y,phi,'linewidth',2) %y,phi_i,y,phi_e,
axis([.4 .6 -.1 .1])
xlabel('Length (cm)','fontsize',16)
ylabel('Potential','fontsize',16)

%test=[y' phi'];
%save('test.dat','-ascii','test');

figure(2)
plot(y,phi_i,y,phi_e)


