
%  -------------------------------------------------------------------
%
%  % find the speed of solitary pulses for  piecewise linear FHN system
%
%   For Chapter 6, Figure 1.8, 1.9 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 

 clear all; close all; clc;
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);  
 

global lam c eps pr x1
% parameters


epslist = [0.5,0.1,0.01];

for jk = 1:3
    eps = epslist(jk);
    clist = [sqrt(eps):0.002:2];
for jj = 1:length(clist)
    c=clist(jj);
% find roots of the polynomial
 
lam = findlam(c, eps);
 pr = pprime(lam);
 
x=[0:.01:10];
 
%find the root by bisection
 % guess a lower bound for the root
a = 0.0001;

% guess an upper bound for the root
b = 10;
   xans(jj)  = bisect(@h,a,b);
   alp(jj) = (1-exp(lam(1)*xans(jj)))/( pr(1));

end
figure(1)
plot(alp,clist)
hold on 
end
figure(1)
clist=[0:.01:2];
alp0 = clist./sqrt(clist.^2+4);
plot((1-alp0)/2,clist,'--')
xlabel('\alpha')
ylabel('c')
box off
legend('boxoff')
legend('\epsilon=0.5','\epsilon=0.1','\epsilon=0.01')

% now plot profiles of the traveling wave solution
% specify c and epsilon
 clist = [0.34,2.66];
 xend=[1,6];
 xlow=[-1.5,-15];
  for j =1:2
      c=clist(j);
eps = 0.1;
lam = findlam(c,eps);
pr = pprime(lam);
 lam(1)*lam(2)*lam(3)*c*eps^2
%find the root by bisection
 
   x1  = bisect(@h,a,b);
   
   alf = (1-exp(lam(1)*x1))/(pr(1))
x = [xlow(j):.01:xend(j)];
[v,w] = profiles(x);
figure(1+j)
plot(x,w,'--',x,v,x,zeros(length(x),1),'k--')
box off
xlabel('x')
legend('boxoff')
legend('v','w','location','northwest')
 formatSpecF = '%5.3f\n';
 
    title(strcat('\alpha = ',sprintf(formatSpecF,alf),', c = ',sprintf(formatSpecF,c)))
  end



%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = h(x)
global lam pr
out = 2- exp(lam(1)*x)+pr(1)*exp(-lam(2)*x)/pr(2) +pr(1)*exp(-lam(3)*x)/pr(3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ppr =pprime(x)
global c eps

ppr = 3*eps^2*x.^2 +2*c*eps*x -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,w] = profiles(x)
global lam pr x1 c
l1 = lam(1);
l2 = lam(2);
l3 = lam(3);

A = -l2*l3*(exp(l1*x1) - 1)/(exp(l1*x1)*(l1 - l3)*(l1 - l2));
    B1 = -l2*l3/((l1 - l3)*(l1 - l2));
    B2 = l3*l1/((l2 - l3)*exp(l2*x1)*(l1 - l2));
    B3 = -l1*l2/((l1*l2 - l1*l3 - l2*l3 + l3^2)*exp(l3*x1));
    C2 = -l1*l3*(-1 + exp(l2*x1))/(exp(l2*x1)*(l2 - l3)*(l1 - l2));
    C3 = l1*l2*(exp(l3*x1) - 1)/((l1*l2 - l1*l3 - l2*l3 + l3^2)*exp(l3*x1))
 
% create w
w = (x>=x1).*A.*exp(l1*x)+ ((x<x1)&(x>=0)).*(1+B1.*exp(l1*x)+B2.*exp(l2*x)+B3.*exp(l3*x)) ...
    +(x<0).*(C2.*exp(l2*x)+C3.*exp(l3*x));  
   
% create v = -c *w'
v = -c*((x>=x1).*A.*exp(l1*x)*l1+ ((x<x1)&(x>=0)).*( l1.*B1.*exp(l1*x)+B2.*l2.*exp(l2*x)+B3.*l3.*exp(l3*x)) ...
    +(x<0).*(C2.*l2.*exp(l2*x)+C3.*l3.*exp(l3*x))); 

end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function lam = findlam(c,eps)

p = [eps^2 eps*c -1 1/c];

r = roots(p);

n1=find(r<0);

lam(1) = r(n1);
k = 1;
for j = 1:3
    if (j ~= n1)
        k = k+1;
        lam(k)=r(j);
    end
end

 end

%%%%%%%%%%%%%%%%%%%%%%%
% this is the bisection algorithm
function root = bisect(feval,a,b)

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
end

