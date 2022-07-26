function main
clear all; close all; clc;

global input Vp Vi Vg E tp ti td Rm a1 C1 C2 C3 C4 C5 Ub U0 Um Rg alpha beta

input = 100;
td=12;
Vp=3; Vi=11; Vg=10; E=0.2; tp=6; ti=100;  Rm=209; a1=6.6; C1=300;
C2=144; C3=100; C4=80; C5=25.86; Ub=72; U0=4; Um=94; Rg=180; alpha=7.5; beta=1.77;

initial = [1 1 1 1 1 ];
tspan = linspace(0,2000,1000);

[tt,sol]=ode45(@sturisfun,tspan,initial);

figure(1)
plot(tt,sol(:,1)/Vp)
figure(2)
plot(tt,sol(:,3)/(10*Vg))  % Output per dL, as in Sturis

%save sturis.dat sol -ASCII

%%
function dydt=sturisfun(t,y)
global Vp Vi E tp ti td input

I=y(1); G=y(2); h1=y(3); h2=y(4); h3=y(5); 

dydt = [eff1(G) - I/tp
        eff5(h3) + input - eff2(G) - eff3(G)*eff4(I)
        (I - h1)/td
        (h1-h2)/td
        (h2-h3)/td
        ];

%%
function out=eff1(z)
global  Vg  Rm a1 C1 
out=Rm/(1+exp(-z/(C1*Vg) + a1));

%%
function out=eff2(z)
global  Vg  C2 Ub 
out=Ub*(1-exp(-z/(C2*Vg)));

%%
function out=eff3(z)
global  Vg  C3
out=z/(C3*Vg);

%%
function out=eff4(z)
global Vi E  ti  C4  U0 Um  beta
kappa=(1/Vi + 1/(E*ti))/C4;
out=U0 + (Um-U0)/( 1+(kappa*z)^(-beta) );

%%
function out=eff5(z)
global Vp  C5  Rg alpha 
out=Rg/(1+exp(alpha*(z/(C5*Vp)-1)));