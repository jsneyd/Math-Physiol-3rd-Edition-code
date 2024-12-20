%  -------------------------------------------------------------------
%
%   Velocity of a motor pulling a cargo.
%
%   For Chapter 15, Section 15.10.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------
clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0,...
   'defaultlinemarkersize',18);


dr = 0.1;  %the ratio of diffusion coefficients D2/D1
w = [0.1:.1:10];
% The hard spring limit
ex=exp(w);
v = 2*(dr/(1+dr))*(ex-1)./(ex+1);

% the soft spring limit.
% use bisection to find wl
n = length(w);
wa = zeros(1,n);
wl = wa;
f = 2*(ex-1)./(ex+1);;
fa = f-dr*wl;
wb = w;
wl = wb;
f = wl.^2.*(1-exp(w-wl))./((exp(w)-1).*(exp(-wl)-1)+wl.*(exp(w-wl)-1));
fb = f-dr*wl;

%now do bisection
for j = 1:25
    wl = (wa+wb)/2;
    f = wl.^2.*(1-exp(w-wl))./((exp(w)-1).*(exp(-wl)-1)+wl.*(exp(w-wl)-1));
    ft = f-dr*wl;
    ftest = ft.*fa;
    indx = find(ftest>=0);
    wa(indx) = wl(indx);
    indx = find(ftest<0);
    wb(indx) = wl(indx);

end
v2 = dr*wl;

plot(w,v, w ,v2,'--')
text(5,.16,'Hard spring','fontsize',18)
text(5,.33,'Soft spring','fontsize',18)
xlabel('\omega_0')
ylabel('v\delta/ D_1')
box off

out=[w' v' v2'];
%save('cargo.dat','-ascii','out')
