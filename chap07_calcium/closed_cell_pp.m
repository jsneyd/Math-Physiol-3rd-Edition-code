% simple closed cell oscillations
function closed_cell
close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 

%parameters
p.Vs = 0.9;
p.Vp=0.1;
 
p.Kbar = 2.e-7;
p.kf=1.11;
 
p.Ks = 0.1;
p.Kp = 0.3;
 
p.gm = 5.5;
p.p=0.4;  %ip3 parameter 
p.ct = 20;
p.k1 = 400;
p.k2 = 0.2;
p.k3 = 400;
p.k4 = 0.2;
p.k5 = 20;
p.km1 = 52;
p.km2 = 0.21;
p.km3 = 377.2;
p.km4 = 0.0289;
p.km5 = 1.64;
p.K1=p.km1/p.k1;
p.K5=p.km5/p.k5;
p.K2 = p.km2/p.k2;
p.K3=p.km3/p.k3;
p.K4=p.km4/p.k4;


 init = [0.5 ,0.2];
dt=0.1;
tend=100;

tspan = [0:dt:tend];
[t,sol] = ode15s(@(t,x)coscrhs(t,x,p),tspan,init);

figure(1)
plot(t,sol(:,1),'r',t,sol(:,2),'b')  % time series
xlabel('Time')
legend('c','y')

 

% nullclines
cmx = p.kf*p.ct/(p.kf*p.gm+1);
c = [0.009:0.001:cmx];
ph1  = (p.km4*p.K1*p.K2+p.km2*p.K4*p.p).*c./(p.K4*p.K2*(p.p+p.K1));
ph2 =  (p.km2*p.p+p.km4*p.K3)/(p.K3+p.p);

 y1=ph1./(ph1  +ph2);
 
 ce=p.ct-p.gm*c;
 
 Jserca = p.Vs*(c.^2-p.Kbar*ce.^2)./(p.Ks^2+c.^2);
Po=Jserca./(p.kf*ce-c);
 

tmp=((p.p+p.K1).*(c+p.K5).*Po.^(1/3))./(p.p*c);
y2=1- tmp;
 

 figure(2)  % phase portrait
semilogx(sol(:,1),sol(:,2),c,y1,'--',c,y2,'b--')
xlabel('c')
ylabel('y')
axis([0 cmx 0 1])
text(1,0.9,'dy/dt=0','fontsize',18)
text(1,0.1,'dc/dt=0','fontsize',18)

function out=coscrhs(t,x,p)

c=x(1); % calcium
y=x(2); %IPR inhibition
ce = p.ct-p.gm*c;

Po = (p.p*c*(1-y)./((p.p+p.K1).*(c+p.K5))).^3; %open probability
ph1  = (p.km4*p.K1*p.K2+p.km2*p.K4*p.p).*c./(p.K4*p.K2*(p.p+p.K1));
ph2 =  (p.km2*p.p+p.km4*p.K3)/(p.K3+p.p);

Jipr = p.kf*Po*(ce-c);
Jserca = p.Vs*(c.^2-p.Kbar*ce^2)/(p.Ks^2+c^2);

  
out(1) = Jipr-Jserca;
out(2) = ph1*(1-y) -ph2*y;

out = out';
 