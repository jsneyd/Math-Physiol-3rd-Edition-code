

%  -------------------------------------------------------------------
%
%   Code for the compound bursting model of Wierschem and Bertram, 2004.
%   Original code was provided by Richard Bertram and modified slightly.
%
%   For Chapter 9, Section 9.1.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by Richard Bertram, James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 

close all
clear all
clc

global vn sn vm sm kd kadp gca vca gkca vk gkatp gk cm f alpha kc v1 tauc v2 taun

 
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 

%parameters
 
% Conductances in pS; currents in fA; Ca concentrations in uM; time in ms
 
gkatp=350; f=0.001; v2=185; 
% exercise parameters:
%gkatp=357; f=0.0013; v2=188; 

tauc=1.2e6;
gca=1200; gk=3000;
gkca=300; taun=16;
vca=25; vk=-75;
 cm=5300; kc=0.1;
 kd=0.3; alpha=2.25e-6;
 sm=12; vm=-20;
sn=5.6; vn=-16;
v1=10; kadp=20;

% initial data

v0=-61.21;
n0=0.000311;
c0=0.07;
atp0=1.99;
adp0=0.0502;

u0 = [v0,n0,c0,atp0,adp0];
 
total=600000;
tstep = 1;
tic
%specify the output points
tspan = [0:tstep:total];
% warning: This problem is sensitive to integrator choice!
[T,S] = ode23s(@deRHS,tspan, u0, odeset('maxstep',10));  
toc
figure(1)
plot(T/60000,S(:,1))
ylabel('V')
xlabel('t (s)')

yyaxis('right')
plot( T/60000,S(:,4))
ylabel('ATP')
figure(2)
plot(T/1000,S(:,3))
xlabel('t (s)')
ylabel('Ca^{++}')

figure(3)
plot(T/60000,S(:,4))
xlabel('t (s)')
ylabel('ATP')


%%
function s_prime=deRHS(t,s)
global vn sn vm sm kd kadp gca vca gkca vk gkatp gk cm f alpha kc v1 tauc v2 taun

v=s(1);
n=s(2);
c=s(3);
atp=s(4);
adp=s(5);
 
% # Activation variables
ninf = 1/(1+exp((vn-v)/sn));
minf = 1/(1+exp((vm-v)/sm));
omega = 1/(1+(kd/c));

% ATP handling

phi=atp*(1+kadp*adp)^2;

% Ionic currents
ica = gca*minf*(v-vca);
ikca = gkca*omega*(v-vk);
ikatp = (gkatp/atp)*(v-vk);
ik = gk*n*(v-vk);

% Differential Equations

vp = (-ica - ik - ikatp  - ikca)/cm;
np =  (ninf - n)/taun;
cp = -f*(alpha*ica + kc*c);
atpp=(v1-phi)/tauc;
adpp=(phi-v2*adp)/tauc;

s_prime = [vp;np;cp;atpp;adpp];
end