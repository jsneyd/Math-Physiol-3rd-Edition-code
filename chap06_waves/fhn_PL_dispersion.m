% Code to find the dispersion curve for the piecewise linear FHN

%
%   For Chapter 6, Figure 6.13,   of
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
 
global lamm c eps xi2 pr alp

% parameters
 
epslist = [0.1 ,0.01];
alp = 0.1; %target value of alpha

for jk = 1: length(epslist)
    eps = epslist(jk)
    clist =  [0.005:0.0075:2.7];
    P=[];
    cl = [];
 
for jj = 1:length(clist)
    c=clist(jj);
% find roots of the polynomial
 
lam = findlam(c, eps);
lamst(jj,:) = lam;
  prc = [3*eps^2 2*c*eps  -1];
 pr = polyval(prc,lam);
 
 %This is necessary because real parts of Lam(2) and lam(3) are positve
% and can give large exponentials
% real part of lam(1) <0
% real part of lam(2,3) >0
 
lamm(1) = lam(1);
lamm(2) = -lam(2);
lamm(3) = -lam(3);


% find xi2 for this value of c and alpha
 % 
 
     newalp= alproot;
 P = [P,newalp];
 
 cl = [cl,clist(jj)*ones(1,length(newalp))];
end
  for jj = 1:length(P)
xi2 = P(jj);
a=1.e-6;
b = xi2/2- 0.01;
 
  xi1(jj) = bisect(@h,a,b);
 end   
 
figure(1)
plot(P./cl,cl,'*')

xlabel('T')
ylabel('c')
 hold on
%axis([0 10 0 3])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure(2)
 plot(P, cl ,'*' )
 xlabel('\xi_2')
 ylabel('c')
 hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure(3)
  plot(P,xi1,'*')
 xlabel('xi_2')
 ylabel('xi_1')
 hold on
 figure(4)
 plot(clist,real(lamst),clist,imag(lamst))

end
 % The singular dispersion curve
   w0=[0:.01:(1-2*alp)/2];
   apw = alp+w0;
   w1=1-2*alp-w0;
   sp=(1-2*apw)./sqrt(apw-apw.^2);
   T = log(((1-w0).*w1)./((1-w1).*w0));

   figure(1)
    plot(T,sp,'k--')
    box off
    legend('boxoff')
    legend('\epsilon=0.1','\epsilon = 0.01','asymptotic limit','location','northwest')
   hold off

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

 %%%%%%%%%%%%%%%%%%%%%%%%%%
function f = h0(x)
 global xi2
 xi2 = x;

r = getrat(x/2);
f = r(1)-r(2);

end

 %%%%%%%%%%%%%%%%%%%%%%%%%%
function f = h(x)
 
r = getrat(x);
f = r(1)-r(2);

end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = getrat(x)
global lamm  xi2 pr

 
xi = [x, xi2,xi2-x];
lamX=lamm'*xi;
E=exp(lamX);  % These are all the exponentials needed;  they are all less than 1;

Fac = [1,E(1,1);E(2,1),1;E(3,1),1];
  
rat = (1-E(:,3))./(pr'.*(1-E(:,2)));  
  
r1 = rat'*Fac;
 r = rat(1)*Fac(1,:)+rat(2)*Fac(2,:) +rat(3)*Fac(3,:);
 
end


 %%%%%%%%%%%%%%%%%%%%%%%%%%
 function afc  = getalp(x)
global alp xi2
% x is a guess for xi2
 
xi2 = x;
a=0.01;
 
b = x/2;
  
  x1 = bisect(@h,a,b);
  

r = getrat(x1);
afc = r(1)+alp;
 
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

% if not, the algorithm fails to find a root.

N = 20;  % number of iterates
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

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function art = alproot(x)
% this part is slow
% there may be multiple roots

 % find the range of xi2
 xi2min = bisect(@h0,0.01,20);

 xi2list=[xi2min:.01:20];

 for jj = 1:length(xi2list)
     xi2=xi2list(jj);
     f(jj)=getalp(xi2);
 end

 tst = f(1:end-1).*f(2:end);
 ndx=find(tst<0);
 if isempty(ndx)
     art=[];
 else
 art=xi2list(ndx);
 end
end

