
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



global lam c eps pr 
% parameters

c = 0.3
eps = 0.05;

% find roots of the polynomial

p = [eps^2 eps*c -1 1/c];

r = roots(p)

n1=find(r<0);

lam(1) = r(n1);
k = 1;
for j = 1:3
    if (j ~= n1)
        k = k+1;
        lam(k)=r(j);
    end
end

lam
pr = pprime(lam)
s = [0.001:.01:.99];
plot(s,h(s))
fsolve(@h,0.5)

%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = h(s)
global lam c eps pr
out = 2-s+pr(1)*exp(-lam(2)*log(s)/lam(1))/pr(2) +pr(1)*exp(-lam(3)*log(s)/lam(1))/pr(3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ppr =pprime(x)
global c eps

ppr = 3*eps^2*x.^2 +2*c*eps*x -1;
end
