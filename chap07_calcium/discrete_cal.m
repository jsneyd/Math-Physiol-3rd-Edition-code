%this file makes some plots for the discrete release standing waves

%two parameters
Abyk =5;
beta = [0.01:.1:5];

astar = 1-sqrt(1/(1+Abyk))
Af = Abyk.*beta.*sinh(beta);
l0 = cosh(beta);
lf = l0+Af/2;
mu0 = l0-sqrt(l0.^2-1);
muf = lf-sqrt(lf.^2-1);

C = Af./(Af+2*l0-2);

a = C.*(muf-1)./(muf-1./mu0);
b = C.*(1-1./mu0)./(muf-1./mu0);
bf = C-b.*muf;
plot(beta,a,'r',beta,bf,'b','linewidth',2)
xlabel('\beta','fontsize',16)
ylabel('c^*/c_e','fontsize',16)
text(2,.4,'standing waves','fontsize',16)
text(.5, .1,'traveling waves','fontsize',16)

output=[beta' a' bf'];

save discrete.dat output -ascii
