
%    -------------------------------------------------------------------
%
%     Solve the Sturis model for pulsatile insulin secretion.
%
%     For Chapter 16, Section 16.7.3 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------


function sturis
clear all; close all; clc;
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0);

global input Vp Vi Vg E tp ti td Rm a1 C1 C2 C3 C4 C5 Ub U0 Um Rg alpha beta

% pick a glucose input level
input = 216;
%input = 108;
inputlist=[108,216];
for ic = 1:2
    input = inputlist(ic);


Vp=3; Vi=11; Vg=10; E=0.2; tp=6; ti=100; td=12; Rm=209; a1=6.6; C1=300;
C2=144; C3=100; C4=80; C5=25.86; Ub=72; U0=4; Um=94; Rg=180; alpha=7.5; beta=1.77;

initial = [1 1 1 1 1 1];
tspan = linspace(0,2000,1000);

[tt,sol]=ode45(@sturisfun,tspan,initial);
formatSpecF = '%6.1f\n';
figure(ic)
plot(tt,sol(:,1)/Vp)
xlabel('t (min)')
ylabel('I_i (mU/ml)')
yyaxis right
plot(tt,sol(:,3)/(10*Vg))  % Output per dL, as in Sturis
%xlabel('t (min)')
ylabel('G (mg/dl)')
%save sturis.dat sol -ASCII
%title('')
box off
 title(strcat('I_G = ',sprintf(formatSpecF,input)),'fontsize',20)
end
%%
function dydt=sturisfun(t,y)
global Vp Vi E tp ti td input

Ip=y(1); Ii=y(2); G=y(3); h1=y(4); h2=y(5); h3=y(6);

dydt = [eff1(G) - (Ip/Vp-Ii/Vi)*E - Ip/tp
        (Ip/Vp-Ii/Vi)*E - Ii/ti
        eff5(h3) + input - eff2(G) - eff3(G)*eff4(Ii)
        (Ip - h1)/td
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
