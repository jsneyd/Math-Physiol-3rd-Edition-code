% code to model postsynaptic respone to presynaptic voltage action potential

function full_model
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 
%parameters
global Vr k10  k20  z1  n ce PCa s0 
global FRT  F  kout  radius 
global   kp1  km1 kp2  km2  kp3  km3 kp4   km4 
global   bigA  bigB   bigN smalla  smallb k1m  k2m  ke 
global Vs  gr  gs0 Cm Istim Tperiod  rho

% FHN parameters
Vr=-70;

% Llinas parameters
 k10=2;
 k20=1;
 z1=1;
 n=5;
 %ci=0.1;
 ce=4.0e4;
 PCa=3.0e-2;
 s0=100;

%Voltage-type parameters
 FRT=.0386;
 F=96490.0 ; 

% Calcium parameters. Radius in microns. Factor converts current to calcium flux.
  kout=1; 
  radius=20;

%%!factor=(3/radius)*(1/(2*F))

% Calcium-stimulated secretion parameters
 kp1=3.75e-3;
 km1=4e-4;
 kp2=2.5e-3; 
 km2=1e-3;
 kp3=5e-4;
 km3=0.1;
 kp4=7.5e-3;
 km4=10;
 rho = 3e6;

% Magleby parameters
 bigA=0.008;
 bigB=1.43;
 bigN=10;
smalla=0.00315; 
smallb=1;
 k1m=1000; 
 k2m=500;
 ke=10;% rate of ACh degradation

% Postsynaptic voltage parameters
 Vs=-15;
 gr=10;
 gs0 = 10;
 Cm=1;


% there are twelve equations
 Istim=3;
 Tperiod = 30;

tt=[0:.01:5];
init=[-50, 0,0,0,0,0,0,0,0,0,0,-70];
[T,sol] = ode15s(@(t,S)rhs(t,S),tt,init);
v1=sol(:,1);
v2=sol(:,12);
op=sol(:,3);
c = sol(:,4);
a = sol(:,11);
x = sol(:,9);
y=sol(:,10);  

figure(1)

 plot(T,v1)
 legend('boxoff')
 legend('V1','V2')
xlabel('t (ms)')

figure(2)
plot(T,c)
legend('c','op')
xlabel('t (ms)')

figure(3)
plot( T,op)
legend('c','op')
xlabel('t (ms)')

figure(4)
plot(T,a,T,x,T,y)
legend('a','x','y')
xlabel('t (ms)')


%%
function out = rhs(t,sol)
global Vr k10  k20  z1  n ce PCa s0 
global FRT  F  kout  radius 
global   kp1  km1 kp2  km2  kp3  km3 kp4   km4 
global   bigA  bigB   bigN smalla  smallb k1m  k2m  ke 
global Vs  gr  Cm gs0 rho

%there are 12 variables
v1=sol(1); %presynaptic potential
w = sol(2);
o = sol(3); % open probability
c = sol(4); % calcium concentration
o1 = sol(5);
o2 = sol(6);
o3 = sol(7);
o4 = sol(8);
x = sol(9);
y = sol(10);
a = sol(11);
v2 = sol(12);
R=o1*o2*o3*o4;

%%%%%%%%%%%%%%%%%%%%%
% Llinas functions  
k1=k10*exp(z1*v1*FRT);
fac=2*v1*FRT;
 
if abs(fac)<1.e-6
    g=1/(1-fac/2);
else
g = fac/(1-exp(-fac));
end

j=PCa*2*F*g*(c-ce*exp(-fac));

ICa=j*s0*((o)^n);

% calcium-stimulated secretion functions
tau1=1/(kp1*c+km1);
tau2=1/(kp2*c+km2);
tau3=1/(kp3*c+km3);
tau4=1/(kp4*c+km4);

%Magleby functions
alpha=bigB*exp(bigA*v2);
beta=smallb*exp(smalla*v2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% FHN odes for the presynaptic voltage
v1p=100*(0.0001*(v1-Vr)*(70-(v1-Vr))*((v1-Vr)-7)-w);
wp =0.25*(v1-Vr-5*w);

%Llinas odes
op=k1*(1-o) - k20*o;

%Calcium odes. Conversion factors are included explicitly, for clarity

 kout=1;
 radius=20;
factor=(3/radius)*(1/(2*F));

cp= -ICa*factor - kout*c;

%Calcium-stimulated secretion odes
o1p=kp1*c - o1/tau1;
o2p=kp2*c - o2/tau2;
o3p=kp3*c - o3/tau3;
o4p=kp4*c - o4/tau4;

% Magleby odes
xp=-alpha*x + beta*y;
yp=alpha*x + k1m*a*(bigN-x-y) - (beta+k2m)*y;
ap=rho*o1*o2*o3*o4 - ke*a - k1m*a*(bigN-x-y) + k2m*y;

% Postsynaptic voltage ode
v2p= (1/Cm)*(-gr*(v2-Vr) - gs0*x*(v2-Vs));

out  =[v1p,wp,op,cp,o1p,o2p,o3p,o4p,xp,yp,ap,v2p]';
 

 
