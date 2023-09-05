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
global FbyRT  F  kout  factor  
global   kp1  km1 kp2  km2  kp3  km3 kp4   km4 
global   bigA  bigB   bigN smalla  smallb k1m  k2m  ke 
global Vs  gr  gs0 Cm Istim Tperiod  rho

tmax = 15;  %this is the length of the simulation
% FHN parameters
Vr=-70;  % resting potential 
Tperiod =tmax;  % period of stimulus
 Istim=100; % stimulus amplitude

 % Llinas parameters
 k10=2;
 k20=1;
 z1=1;
 n=5;
 
 ce=4e4; % units micro molar = 40 mM
 PCa=2.0e-2; % mu m ms^{-1}
 s0=20;

%Voltage-type parameters
F=96490; % Faraday's constant
R=8.315;  % Gas constant
Tp=300; % Temperature in degrees Kelvin
 FbyRT = F/(R*Tp*1000); % units of mV^(-1)
  
% Calcium parameters. Radius in microns. Factor converts current to calcium flux.
  
factor=(3/20)*(1/(2*F));
kout = 3; % calcium withdrawal rate

% Calcium-stimulated secretion parameters for o_j dynamics
 kp1=3.75e-3;
 km1=4e-4;   
 %km1=1e-4;%these rates are extremely slow.   
 kp2=2.5e-3; 
 km2=1e-3;   
 %km2 = 5e-2; %these rates are extremely slow.  
 kp3=5e-4;
 km3=0.1;
 kp4=7.5e-3;
 km4=10;
 rho = 3e6;

% Magleby parameters
 bigA=0.008;
 bigB=1.43;
 bigN=10;  %dimensionless
smalla=0.00315; 
smallb=1;
 k1m=1000; %ms^(-1) 
 k2m=500; %mM ms^(-1)
 ke=10;% rate of ACh degradation

% Postsynaptic voltage parameters
 Vs=-15; 
 gr=10;
 gs0 = 10;
 Cm=1;


% there are twelve equations

tt=[0:.01:tmax];
init=[-70, 0,0.11,1.4e-6,0.14,0.05,0.0015,0,0,0,0,-70];
[T,sol] = ode15s(@(t,S)rhs(t,S),tt,init);


v1=sol(:,1);
w = sol(:,2);
op=sol(:,3);
c = sol(:,4);
o1=sol(:,5);
o2=sol(:,6);
o3=sol(:,7);
o4=sol(:,8);
Pr = o1.*o2.*o3.*o4;

x = sol(:,9);

y=sol(:,10); 
a = sol(:,11);
v2=sol(:,12);

 figure(1)
 plot(T,v1)
xlabel('t (ms)')
ylabel('V_1 (ms)')
 
figure(2)
plot( T,op)
ylabel('o')
xlabel('t (ms)')

figure(3)
plot(T,c)
 ylabel('c (\muM)')
xlabel('t (ms)')
 
figure(4)
plot(T,o1,T,o2,T,o3,T,o4 )
legend('boxoff')
legend('o_1','o_2','o_3','o_4')
xlabel('t (ms)')

figure(5)
plot(T,Pr)
ylabel('P_R') 
xlabel('t (ms)')
 
figure(6)
 plot(T,a )
 xlabel('t (ms)')
ylabel('a (mM)')

figure(7)
 plot(T,x,T,y,'--')
 legend('boxoff')
 legend('x','y')
 xlabel('t (mV)')
 
 figure(8)
 plot(T,v1,T,v2,'--')
 legend('boxoff')
 legend('presynaptic,','postsynaptic')   
xlabel('t (mV)')
ylabel('Membrane Potential (ms)')

%%
function out = rhs(t,sol)
global Vr k10  k20  z1  n ce PCa s0 
global FbyRT  F  kout  factor
global   kp1  km1 kp2  km2  kp3  km3 kp4   km4 
global   bigA  bigB   bigN smalla  smallb k1m  k2m  ke 
global Vs  gr  Cm gs0 rho Tperiod Istim

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
k1=k10*exp(z1*v1*FbyRT);
fac=2*v1*FbyRT;
 
if abs(fac)<1.e-6  % in case fac is close to zero:
    g=1+fac/2+fac^2/12-fac^4/720;

else
g = fac/(1-exp(-fac));
end

jj=PCa*2*F*g*(c-ce*exp(-fac));

ICa=jj*s0*((o)^n);

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
 
s=t-fix(t/Tperiod)*Tperiod -1;

Iapp = Istim*1./cosh(5*s);

v1p=100*(0.0001*(v1-Vr)*(70-(v1-Vr))*((v1-Vr)-7)-w) +Iapp;
wp =0.25*(v1-Vr-5*w);

%Llinas odes
op=k1*(1-o) - k20*o;

%Calcium odes. Conversion factors are included explicitly, for clarity
 
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
 

 
