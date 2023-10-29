%determine stability for the FHN cable equation with pl f(v)
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0);
%parameters
gam = 0.1;
alp  = 0.2;

lam = sqrt(1+1/gam);
mu = sqrt(1/gam-1);

v0=alp/(1+gam);
v1=(alp-gam/2)/(1-gam);
%region 1
Imx = lam*(1/4-v0);

I=[0:.01:Imx];
vmx = v0+I/lam;

figure(1)
plot(I,vmx)
hold on


%region 2

A=(-4*v1 + 1)*mu  - 4*(v0 - 1/4)*lam;
B= (8*v1 - 6)*mu;
C=4*(v0 - 1/4)*lam +(1- 4*v1)*mu;

D=sqrt(B^2-4*A*C);
Elim=(-B+D)/(2*A);
xlim=log(Elim)/mu;

x1=[0:.001:xlim]; %a guess
 
E=exp(mu*x1);
I = ((-4*v0 + 1)*lam + (-4*v1 + 1)*mu).*E/8 + ((-4*v0 + 1)*lam + mu*(4*v1 - 1))./(8*E);
b = -(4*I./E + 4*v1*mu - mu)./(4*(E + 1./E)*mu); 
c = (4*I.*E - 4*v1*mu + mu)./(4*(E + 1./E)*mu);
V0 = v1+b+c;

 
P=((-4*v1 + 1)*mu  - 4*(v0 - 1/4)*lam)*E.^2 + (8*v1 - 6)*mu*E   + 4*(v0 - 1/4)*lam +(1- 4*v1)*mu;

figure(2)
plot(I,x1)
xlabel('I')
ylabel('x_1')
 hold on

figure(1)
plot(I,V0)
xlabel('I')
ylabel('v(0)')

%region 3
x1=[xlim:0.01:1.5];
E1=exp(mu*x1);
fac=(-4*v1*mu - sqrt(16*v0^2*lam^2 - 8*v0*lam^2 - 16*v1*mu^2 + lam^2 + 8*mu^2) + 3*mu)/(4*v0*lam - 4*v1*mu - lam + mu)
E2=fac*E1;
x2=log(E2)/mu;
 figure(3)
  plot(x1,x2)
 E3=exp(lam*x2);
 v2=(gam+alp)/(gam+1);

  
 % Find I
 v2fac=-4*lam*( E1.^2./E2 - E2).*(E3 - 1./E3);
 v1fac = -4*mu*( E1.^2./E2 - 2*E1  + E2).*(E3 + 1./E3);
E1sqfac = 3* (E3*lam + E3*mu - lam./E3  + mu./E3 )./E2;
E1fac= -2 *mu*(E3 + 1./E3);
E10fac=-3*E2.*(E3*lam - E3*mu - lam./E3 -  mu./E3);
denom=8*(E1.^2./E2 - E2);
I=(v2*v2fac+v1*v1fac+E1.^2.*(E1sqfac) +E1.*E1fac+E10fac)./denom;
    
 figure(2)
 plot( I,x1,'--',I,x2,'--')

 hold off


