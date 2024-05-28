
%  -------------------------------------------------------------------
%
%   Code for the three-variable Chay-Keizer model of electrical bursting
%   in pancreatic beta cells.
%
%   For Chapter 9, Section 9.1.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd
%
%  -------------------------------------------------------------------

function ChayKeizer3v
global   lam vca vk gk gca   vprime vstar gkca Kd
global  vl cm   f alpha kca  gl

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);
 % parameter set

 cm=4.5;
 gca=13;
 gk=12;
 gl=0.04;
 gkca=0.09;
 vca=100;
 vk=-75;
 vl=-40;
  vprime=50;
  vstar=30;
 f=0.01;
 alpha=0.01525;
 Kd=1;
 kca=0.04;
 lam=0.5;


total=10000;
tstep = 1;

%specify the output points
tspan = [0:tstep:total];
%initial data
v0=-54.77;
n0=0.00044;
ca0=0.1075;
 u0 = [v0,n0,ca0];

[T,S] = ode23(@deRHS,tspan, u0, odeset('maxstep',1));

figure(1)
plot(S(:,3),S(:,1))
ylabel('V')
xlabel('Ca^{++}')
figure(2)
plot(T/1000,S(:,1))
xlabel('t (s)')
ylabel('V')

figure(3)
plot(T/1000,S(:,3))
xlabel('t (s)')
ylabel('Ca^{++}')

end % of main


%%
function s_prime=deRHS(t,s)
global  cstar gkcabar lam vca vk gk gca captot vbar vprime vstar gkca Kd
global  vl cm vm sm vh sh a  b c f alpha kca  gl

v=s(1);
n = s(2);
ca=s(3);

% gating functions
alphamca    = -0.1*((v+vprime)-25)/(exp(-((v+vprime)-25)/10)-1);
betamca     = 4*exp(-(v+vprime)/18);
alphahca    = 0.07*exp(-(v+vprime)/20);
betahca     = 1/(exp(-((v+vprime)-30)/10)+1);
alphan      = -0.01*((v+vstar)-10)/(exp(-((v+vstar)-10)/10)-1);
betan      = 0.125*exp(-(v+vstar)/80);
minf =alphamca /(alphamca + betamca );
hinf =alphahca/(alphahca + betahca);

% Ionic currents:
ica  = gca*minf^3*hinf*(v-vca);
ik   = gk*n^4*(v-vk);
ikca = gkca*ca/(ca+Kd)*(v-vk);
il   = gl*(v-vl);

% The differential equations:
vp    = -1/cm*(ica+ik+ikca+il) ;
np    = lam*(alphan*(1-n)-betan*n) ;
cap   = f*(-alpha*ica-kca*ca);

s_prime = [vp;np;cap];

end

