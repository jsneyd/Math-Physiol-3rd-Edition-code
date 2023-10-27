%pl schrodinger equation
global x1 x2
mu =[-1:.01:1];
lam =sqrt(1-mu);
eta=sqrt(1+mu);
x1=1
x2=0.75


  %x1 = atan(eta./lam)./lam;
 
figure(1)
plot(mu,f(mu))

%given x1 and x2 find mu by bisection
function out=f(mu);
global x1 x2
lam =sqrt(1-mu);
eta=sqrt(1+mu);
E=exp(eta*x2);
  out=E.*(cos(lam*x2) .*sin(lam*x1).*(eta.^2-lam.^2) ...
      + 2*cos(lam*x2) .*cos(lam*x1).*eta.*lam ...
      + sin(lam*x2) .*cos(lam*x1).*(lam.^2-eta.^2)) ...
      + (sin(lam*x2) .*cos(lam*x1).*(eta.^2+lam.^2) ...
      - cos(lam*x2) .*sin(lam*x1).*(eta.^2+lam.^2) ...
      + 2*sin(lam*x2) .*sin(lam*x1).*eta.*lam)./E;

end
