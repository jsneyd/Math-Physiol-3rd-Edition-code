
% this is to make  \mu-\tau bifurcation curves for the Layton feedback model


clear

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

global k mu gamma 

%parameters
kbar = 0.5;  
cbar = exp(-kbar);
a = 10;

%
kk = 0;
for k = 0.1:.002:2.1
c = exp(-k);
F = 1+tanh(a*(cbar-c));
Fp = -a*sech(a*(cbar-c))^2;

mu = k*F;
gamma = -k*Fp*c;
 
kk = kk + 1;
 
T(kk) =  fsolve(@Feval,1);
U(kk) = mu;
end
 

dum=[U' T'];
save('tgo.dat','-ascii','dum')
figure(1)
plot(U,T);
axis([0 1.2 0 2])
xlabel('\mu','fontsize',20)
ylabel('\tau','fontsize',20)
text(.45,1,'unstable','fontsize',18)
text(.45,.1, 'stable','fontsize',18)

function out = Feval(tau)
global k mu gamma 
w = pi./(tau +k/(2*mu));
crv = w./(2*sin(k*w/(2*mu)))-gamma;
out = crv;
end

