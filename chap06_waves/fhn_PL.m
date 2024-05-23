
%  -------------------------------------------------------------------
%
%  % find the dispersion curve for  piecewise linear FHN system
%
%   For Chapter 6, Figure 1.8 of
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
   'defaultpatchlinewidth', 0.7); clear all; close all; clc;
 

global lam c eps pr 
% parameters

clist = [0.01:.01:2];
epslist = [0.5,0.1,0.01];

for jk = 1:3
    eps = epslist(jk);
for jj = 1:length(clist)
    c=clist(jj);
% find roots of the polynomial

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

 pr = pprime(lam);
 

%find the root by bisection
 % guess a lower blound for the root
a = 1.e-14;

% guess an upper bound for the root
b = 0.99999;
   sans(jj)  = bisect(@h,a,b);
   alp(jj) = (1-sans(jj))/( pr(1));

end

ndx=find((sans>10*a) & (sans <0.9999));
figure(1)
plot(clist(ndx),sans(ndx))
hold on
 figure(2)
 plot(alp(ndx),clist(ndx))
 hold on
end
figure(2)
alp0 = clist./sqrt(clist.^2+4);
plot((1-alp0)/2,clist,'--')
xlabel('\alpha')
ylabel('c')
box off
legend('boxoff')
legend('\epsilon=0.5','\epsilon=0.1','\epsilon=0.01')


%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = h(s)
global lam   pr
out = 2-s+pr(1)*exp(-lam(2)*log(s)/lam(1))/pr(2) +pr(1)*exp(-lam(3)*log(s)/lam(1))/pr(3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ppr =pprime(x)
global c eps

ppr = 3*eps^2*x.^2 +2*c*eps*x -1;
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

