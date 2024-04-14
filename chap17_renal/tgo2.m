
%    -------------------------------------------------------------------
%
%     Make \mu-\tau bifurcation curves in the Layton model of
%     tubuloglomerular oscillations.
%
%     For Chapter 17, Section 17.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function tgo2

clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

global k mu gamma

%parameters
kbar = 0.5;
cbar = exp(-kbar);
a = 10;
kk = 0;
for k = 0.1:.002:2.1
    c = exp(-k);
    F = 1+tanh(a*(cbar-c));
    Fp = -a*(1/cosh(a*(cbar-c)))^2;

    mu = k*F;
    gamma = -k*Fp*c;
    kk = kk + 1;
    options = optimset('Display','off');
    T(kk) =  fsolve(@Feval,1,options);
    U(kk) = mu;
end

dum=[U' T'];
%save('tgo.dat','-ascii','dum')
figure(1)
    plot(U,T);
    axis([0 1.2 0 2])
    xlabel('\mu','fontsize',20)
    ylabel('\tau','fontsize',20)
    text(.45,1,'unstable','fontsize',18)
    text(.45,.1, 'stable','fontsize',18)


end % of main


%%
function out = Feval(tau)
    global k mu gamma
    w = pi./(tau +k/(2*mu));
    crv = w./(2*sin(k*w/(2*mu)))-gamma;
    out = crv;
end

