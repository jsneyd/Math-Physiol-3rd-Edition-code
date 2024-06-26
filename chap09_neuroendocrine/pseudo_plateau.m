% ------------------------------------
% Pseudo-plateau bursting. A modification of the Chay-Keizer model.
% Modified from the original code of Teka, Bertram, etc.

% Used for Keener and Sneyd, Mathematical Physiology, Chapter 9, Section 9.4.

% Variables:
%    V -- membrane potential
%    n -- delayed rectifier activation variable
%    c -- cytosolic calcium concentration
% ------------------------------------

function pseudo_plateau

clear all
close all
clc

global  gkatp  vk vca cm gk gca  vm  sm  vn sn taun  gkca  kd  kpmca f  alpha
  set(0,   ...
        'defaultaxesfontsize', 20,   ...
        'defaultaxeslinewidth', 2.0, ...
        'defaultlinelinewidth', 2.0);

% parameters

gkatp=180;
vk=-75;
vca=25;
cm=5300;
gk=2700;
gca=1000;
vm=-20;
sm=12;
vn=-12;
sn=5;
taun=18.7;
gkca=400;
kd=0.3;
kpmca=0.18;
f=0.01;
alpha=4.50e-6;

% Initial conditions
v0=-65 ;
n0=0 ;
c0=0.1 ;


init=[v0,n0,c0];

total=10000;
tstep = 0.01;
%specify the output points
tspan = [0:tstep:total];
[T,S] = ode23s(@deRHS,tspan, init, odeset('maxstep',10));

figure(1)
plot(T /1000,S(:,1))
ylabel('V (mV)')
xlabel('t (s)')
box off
figure(2)
plot(S(:,3),S(:,1))
ylabel('V (mV)')
xlabel('Ca^{++}')
box off 
end % of main


%% ode equations
function s_prime=deRHS(t,sol)
global  gkatp  vk vca cm gk gca  vm  sm  vn sn taun  gkca  kd  kpmca f  alpha

%there are three variables
v=sol(1);
n=sol(2);
c = sol(3);

% steady state functions
ninf = 1/(1+exp((vn-v)/sn));
minf = 1/(1+exp((vm-v)/sm));

% Ikca
Ikca = gkca/(1+(kd/c)^3)*(v-vk);

% Calcium Handling
% ICa
Ica = gca*minf*(v-vca);

% Ik
Ik = gk*n*(v-vk);

% Ikatp
Ikatp = gkatp*(v-vk);

% Ca fluxes
Jmem = -(alpha*Ica + kpmca*c);

% equations
dv=-(Ik + Ica + Ikca + Ikatp)/cm;
dn=(ninf-n)/taun;
dc = f*Jmem;

s_prime=[dv;dn;dc];

end

