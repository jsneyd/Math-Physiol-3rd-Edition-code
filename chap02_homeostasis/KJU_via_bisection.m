function KJUtest  % solve the KJU model using bisection

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
  
global  Ks Cs z rhon rhok Nm

% set parameters
 Ks = 2.5;
 Cs = 122.5;
  z = -2;
  rhon = 1.501;
  rhok = 1.01;
  RTbyF = 25.8;
Nmset = [10:5:200];

for j = 1:length(Nmset)
    Nm = Nmset(j);
     V(j) = bisect_function;
end

mu= KJUans(V);
 
figure(2)
plot(Nmset,mu,'linewidth',2)
xlabel('mucosal [Na^+], N_m (mM)')
ylabel('\mu (scaled cell volume)')

figure(1)
plot(Nmset,RTbyF*V,'linewidth',2)
xlabel('mucosal [Na^+], N_m (mM)')
ylabel('V(mV)')

function root=bisect_function
 % a bisection root finder algorithm
% 
% guess a lower bound for the root
a = -3.5;

% guess an upper bound for the speed
b = -0.01;

%call the function bisect(a,b)

root = bisect(a,b);
 
function root = bisect(a,b)

ul = a;
fl = feval(ul);
uu = b;
fu = feval(uu);

% we make the assumption, without checking, that 
% fu*fl<0

% if not, the algoritm fails to find a root.

N = 25;  % number of iterates
% the main bisection algorithm
for j = 1:N
u = (ul+uu)/2;
fc = feval(u);
ftest = (fc*fl>=0);
ul = ftest*u+(1-ftest)*ul;
fl = ftest*fc + (1-ftest)*fl;

uu = (1-ftest)*u+ ftest*uu;
fu = (1-ftest)*fc + ftest*fu;
end
root = u;

function f=feval(v)
global Nm Ks Cs z rhon rhok

%parameters
 Ns = Cs - Ks;
 dum = exp(-v);
 Ci = Cs./dum;

Ni= v.*Nm.*dum./(v + 3*rhon.*(1-dum)) ;
Ki =  Ks*dum + 2*rhok*Ni.*(1-dum)./v;
mu = -z./(Ni+Ki-Ci) ;
f=Ni+Ki+Ci + 1./mu - (Ns+Ks+Cs);

function mu = KJUans(v)
global Nm Ks Cs z rhon rhok
 Ns = Cs - Ks;
 dum = exp(-v);
 Ci = Cs./dum;

Ni= v.*Nm.*dum./(v + 3*rhon.*(1-dum)) ;%Ni
Ki =  Ks*dum + 2*rhok*Ni.*(1-dum)./v; %Ki
mu = -z./(Ni+Ki-Ci) ; %mu
