
% this is to produce \mu-\tau bifurcation curves for the Layton feedback model


clear

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
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
% use bisection to find tau

tl = 0.01;

tau = tl;
w = pi./(tau +k/(2*mu));
crv = w./(2*sin(k*w/(2*mu))) -gamma;
Fl = crv;

tr = 10;
tau = tr;
w = pi./(tau +k/(2*mu));
crv = w./(2*sin(k*w/(2*mu))) -gamma;
Fr = crv;

if(Fr*Fl<0)

for j = 1:30
    tau = (tr+tl)/2;
    
w = pi./(tau +k/(2*mu));
crv = w./(2*sin(k*w/(2*mu)))-gamma;

if(crv*Fl)>0
tl = tau;
Fl = crv;
else
tr = tau;
Fr = crv;
end
end
kk = kk + 1;
T(kk) = tau;
U(kk) = mu;
end
end

dum=[U' T'];
save('tgo.dat','-ascii','dum')

plot(U,T);
axis([0 1.2 0 2])
xlabel('\mu','fontsize',16)
ylabel('\tau','fontsize',16)
text(.45,1,'unstable','fontsize',16)
text(.45,.1, 'stable','fontsize',16)


