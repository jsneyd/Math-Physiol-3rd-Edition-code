% the beginning of code to find dispersion curve for piecewise linear FHN

%
%   For Chapter 6, Figure 1.12,   of
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
 

global lam c eps  

% parameters


epslist = [ 0.1,0.01];

for jk = 1:length(epslist)
    eps = epslist(jk);
    clist = [sqrt(eps):0.002:2];
for jj = 1:length(clist)
    c=clist(jj);
% find roots of the polynomial
 
lam = findlam(c, eps);
 
end
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


 %%%%%%%%%%%%%%%%%%%%%%%%%%
function out = h(x)
global lam alp c 
x1=x(1);
x2=x(2);

%two equations to solve
l1=lam(1);
l2=lam(2)
l3=lam(3);

out(1)=-c*(-(exp(l1*x2) - exp(l1*x1))*l2*l3*l1/((exp(l1*x2)*l1 - exp(l1*x2)*l3 - l1 + l3)*exp(l1*x1)*(l1 - l2)) ...
    + l1*l3*(exp(l2*x2) - exp(l2*x1))*l2/((exp(l2*x2)*l2 - exp(l2*x2)*l3 - l2 + l3)*exp(l2*x1)*(l1 - l2)) ...
    - l1*l2*(exp(l3*x2) - exp(l3*x1))*l3/((l1*l2 - l1*l3 - l2*l3 + l3^2)*(exp(l3*x2) - 1)*exp(l3*x1))) - alp;
 

out(2) = -c*(-(exp(l1*x2) - exp(l1*x1))*l2*l3*l1/((exp(l1*x2)*l1 - exp(l1*x2)*l3 - l1 + l3)*(l1 - l2))  ...
    + l1*l3*(exp(l2*x2) - exp(l2*x1))*l2/((exp(l2*x2)*l2 - exp(l2*x2)*l3 - l2 + l3)*(l1 - l2)) ...
    - l1*l2*(exp(l3*x2) - exp(l3*x1))*l3/((l1*l2 - l1*l3 - l2*l3 + l3^2)*(exp(l3*x2) - 1))) - alp ;
 
end