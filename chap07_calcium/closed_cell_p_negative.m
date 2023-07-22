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

p.Vplc = 0.005;
p.ct = 2;

 
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
%p.K5=p.km5/p.k5;
p.K5=0.2;
p.K2 = p.km2/p.k2;
p.K3=p.km3/p.k3;
p.K4=p.km4/p.k4;

p.gm = 5.5;
p.n = 2;
p.Kbar = 0;
p.kf=7.4;
p.Ki=0.3;
p.kf2=0.00148;
 
p.k3K=0.1;
p.K3K = 0.4;
p.k5P=0;
p.Kplc=0;
p.Vserca = 0.25; 
p.Kserca = 0.1;


p.tauy=6.6;

 init = [0.05 ,0.2,0.2];
dt=0.05;
tend=200;

tspan = [0:dt:tend];
[t,sol] = ode15s(@(t,x)coscrhs(t,x,p),tspan,init);

figure(1)
plot(t,sol(:,1),'r',t,sol(:,2),'b',t,sol(:,3))  % time series
xlabel('Time')
legend('c','y','p')

  
function out=coscrhs(t,x,p)

c=x(1); % calcium
y=x(2); %IPR inhibition
P = x(3); %IP3

ce = p.gm*(p.ct-c);

Po = (P*c*(1-y)./((P+p.K1).*(c+p.K5))).^3; %open probability

Jipr = (p.kf*Po+p.kf2)*(ce-c);
Jserca = p.Vserca*(c.^2-p.Kbar*ce^2)/(p.Kserca^2+c^2);
 
out(1) = Jipr-Jserca;
out(2) =  (-1+(1-y)*(p.Ki+c)/p.Ki)/p.tauy;
out(3) =  p.Vplc*c^p.n/(p.Kplc^2+c^p.n)-(p.k3K+p.k5P)*c.^2/(p.K3K^2+c^2)*P;

out = out';
 