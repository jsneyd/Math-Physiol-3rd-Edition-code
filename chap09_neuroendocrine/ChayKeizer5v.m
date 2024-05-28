%  -------------------------------------------------------------------
%
%   Code for the five-variable Chay-Keizer model of electrical bursting
%   in pancreatic beta cells.
%
%   For Chapter 9, Section 9.1.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd
%
%  -------------------------------------------------------------------

function chayKeizer5v
global vprime vstar gca vca gk gkca Kd vk gl vl cm k1 kca f

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);


% Parameters:
% conductances in mS/cm^2, ca in uM, v in mV, t in msec
  cm=4.5; gca=13; gk=12; gl=0.04; gkca=0.09;
  vca=100; vk=-75; vl=-40;
 vprime=50; vstar=30;
 f=0.01; k1=0.01525;
  Kd=1; kca=0.04;

  total=15000;
tstep = 1;

%specify the output points
tspan = [0:tstep:total];
%initial data

v0    = -54.774
ca0   = 0.10749
mca0  = 0.027532
hca0  = 0.086321
n0    = 0.00044035

 u0 = [v0,mca0,hca0,n0,ca0];


[T,S] = ode15s(@deRHS,tspan, u0, odeset('maxstep',1));

figure(3)
plot(S(:,5),S(:,1))
ylabel('V')
xlabel('Ca^{++}')
figure(1)
plot(T/1000,S(:,1))
xlabel('t (s)')
ylabel('V')

figure(2)
plot(T/1000,S(:,5))
xlabel('t (s)')
ylabel('Ca^{++}')

%%
function s_prime=deRHS(t,s)
global vprime vstar gca vca gk gkca Kd vk gl vl cm k1 kca f
v=s(1);
mca=s(2);
hca=s(3);
n=s(4);
ca=s(5);

% auxiliary variables
alphamca    = -0.1*((v+vprime)-25)/(exp(-((v+vprime)-25)/10)-1);
betamca    = 4*exp(-(v+vprime)/18);
alphahca    = 0.07*exp(-(v+vprime)/20);
betahca     = 1/(exp(-((v+vprime)-30)/10)+1);
alphan      = -0.01*((v+vstar)-10)/(exp(-((v+vstar)-10)/10)-1);
betan       = 0.125*exp(-(v+vstar)/80);

%Ionic currents:
ica  = gca*mca^3*hca*(v-vca);
ik   = gk*n^4*(v-vk);
ikca = gkca*ca/(ca+Kd)*(v-vk);
il   = gl*(v-vl);

% The differential equations:
vp    = -1/cm*(ica+ik+ikca+il) ;
mcap  = (alphamca *(1-mca)-betamca *mca);
hcap  = (alphahca *(1-hca)-betahca *hca) ;
np    = (alphan *(1-n)-betan *n) ;
cap   = f*(-k1*ica-kca*ca);

s_prime = [vp;mcap;hcap;np;cap];
