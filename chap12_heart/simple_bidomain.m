L = 10;
x = [0:.1:L];

rc = 0.5;
lambda = 1;
I = 1;

A = rc*I/lambda;
C  = cosh(lambda*L);
S = sinh(lambda*L);

B = rc*I*(1-C)/(lambda*S);

phi = A*sinh(lambda*x) + B*cosh(lambda*x);
phi_i = -rc*I*x+phi;

figure(1)

plot(x,phi,x,phi_i);


