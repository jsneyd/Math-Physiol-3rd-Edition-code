%the 3 variable Chay Keizer model for beta cells
function CheyKeizer
global  cstar gkcabar lambda vca vk gk gca captot vbar 
global  vn sn vm sm vh sh a  b c f alpha kca kd ibar treset width auto

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 
 cstar=0.7;
 gkcabar=30000;
 lambda=1.6;
 vca=110;
 vk=-75;
gk=2500;
gca=1400;
captot=5309.3;
vbar=-75;
vn=-15;
 sn=5.6;
 vm=4;
 sm=14;
 vh=-10;
 sh=10;
 a=65;
 b=20;
 c=60;
f=0.001;
alpha=4.506e-6;
kca=0.03;
kd=100.0;
ibar=0;
treset=10000;
width=1000;

ibar=0;
treset=10000;
width=1000;
auto=0;


total=15000;
tstep = 1;

%specify the output points
tspan = [0:tstep:total];
%initial data
v0=-62;
n0=0.001;
ca0=0.55;
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
ylabel('V')

function s_prime=deRHS(t,s)
global  cstar gkcabar lambda vca vk gk gca captot vbar 
global  vn sn vm sm vh sh a  b c f alpha kca kd ibar treset width auto

v=s(1);
n = s(2);
ca=s(3);

phik  = 1/(1+exp((vn-v)/sn));
phica  = 1/(1+exp((vm-v)/sm));
phih  = 1/(1+exp((v-vh)/sh));
ica  = gca*phica *phih *(v-vca);
taun  = c/(exp((v-vbar)/a)+exp((vbar-v)/b));
iapp  = ibar*( (t - treset>0) -  (t - treset - width>0));
gkca = gkcabar*ca/(kd + ca);


vp =  (-ica  + gk*n*(vk-v) + gkca*(vk-v) + iapp )/captot;
np =  lambda*(phik -n)/taun ;
cap = auto*(cstar - ca) + (1-auto)*(f*(-alpha*ica  - kca*ca));

s_prime = [vp;np;cap];
% % @ , meth=cvode, toler=1e-6, atoler=1e-6, dt=5
% % @ xhi=1, ylo=-70, yhi=-10, maxstor=20000, bound=100000
% % @ xplot=ca yplot=v
% % aux curpA = (-ica(v) + gk*n*(vk-v) + gkca*(vk-v))/1000
% % aux icaout = -ica(v)/1000
% % aux ikout = gk*n*(vk-v)/1000
% % aux iapp = iapp(t)
% % 
