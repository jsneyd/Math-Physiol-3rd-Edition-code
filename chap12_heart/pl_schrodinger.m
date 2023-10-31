%pl schrodinger equation
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 

 
global x1 x2
mu =[-1:.01:1];  %mu is the proposed eigenvalue
 

lam =sqrt(1-mu);
eta=sqrt(1+mu);
  x1 = atan(eta./lam)./lam;
 
figure(1)
plot(x1,mu)

  stop
 x1=1
x2=0.75

% guess a lower bound
a = -1;

% guess an upper bound
b = 0.99;

%this is the answer - call the function bisect(a,b,@fval)

r = bisect(a,b,@fval)

figure(1)
plot(mu,fval(mu), r,fval(r),'*')


function root = bisect(a,b,feval)

ul = a;
fl = feval(ul);
uu = b;
fu = feval(uu);

% we make the assumption, without checking, that 
% fu*fl<0

% if not, the algoritm fails to find a root.

N = 20;  % number of iterates
% the main bisection algorithm
for j = 1:N
uc = (ul+uu)/2;
fc = feval(uc);
ftest = (fc.*fl>0);
ul = ftest.*uc+(1-ftest).*ul;
fl = ftest.*fc + (1-ftest).*fl;

uu = (1-ftest).*uc+ ftest.*uu;
fu = (1-ftest).*fc + ftest.*fu;
end
root = uc;

end
%given x1 and x2 find mu by bisection
function f=fval(mu);
global x1 x2
lam =sqrt(1-mu);
eta=sqrt(1+mu);
E=exp(eta*x2);
  f=E.*(cos(lam*x2) .*sin(lam*x1).*(eta.^2-lam.^2) ...
      + 2*cos(lam*x2) .*cos(lam*x1).*eta.*lam ...
      + sin(lam*x2) .*cos(lam*x1).*(lam.^2-eta.^2)) ...
      + (sin(lam*x2) .*cos(lam*x1).*(eta.^2+lam.^2) ...
      - cos(lam*x2) .*sin(lam*x1).*(eta.^2+lam.^2) ...
      + 2*sin(lam*x2) .*sin(lam*x1).*eta.*lam)./E;

end
