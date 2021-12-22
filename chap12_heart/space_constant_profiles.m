%this  file calculates the potentials in a passive cardiac cable.

%The sawtooth potential

L = 0.012;
qi = 5.47e-3;
qe = 0.5 *qi;
Q = qi+qe;
sQ = sqrt(Q);
lambda = sQ/L;
lambda_g = 0.09;
mu = exp(-L/lambda_g);

E = exp(sQ);

x = [0:.05:1]*L;

phi = (mu-1/E)*exp(lambda*x) + (mu - E)*exp(-lambda*x);

Pk = 2*qe*(mu-1/E)*(mu-E)/(mu-1);
phi_i = phi*qi+Pk;
phi_e = -phi*qe + Pk;
Vi = [];
Ve = [];
y = [];

sc =1.e3;
for j = 0:15;
    
    Vi = [Vi,sc*mu^j*phi_i];
    Ve = [Ve,sc*mu^j*phi_e];
    y = [y,j*L+x];
end
plot(y,Ve,y,Vi,y,Vi-Ve,'linewidth',2);

xlabel('Length(cm)', 'fontsize',16)
ylabel('Potential', 'fontsize',16)
text(.05,0,'V_e','fontsize',16)
text(.05,-.9,'V_i','fontsize',16)
text(.05,-1.4,'V','fontsize',16)

test=[y' Vi' Ve'];
save('test.dat','-ascii','test')



