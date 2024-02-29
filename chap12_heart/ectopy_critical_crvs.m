%  -------------------------------------------------------------------
%
%   Ectopic focus oscillations Hopf curves
%
%   For Chapter 12, Section 12.4.4 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function main

clear all
close all
clc

set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 1.5, ...
'defaultlinelinewidth', 2.0)

global alps gam Dscal A x

% parameters
gam = 0.1;
b = 0.5;

N=50; % number of interior grid points
L = 5; % length of the domain
dx = L/(N);
x=[1:N]'*dx;
Dl = 1;

Dscal = Dl/dx^2;
for j = 1:15
    alf = j*0.02;
    Alf(j) =alf;

    for k = 1:20
        scale=k*0.125;
        Sc(k) = scale;

        alps=alf +(b-alf)*(exp(-(x/scale).^2));
        A = -2*diag(ones(N,1)) + diag([2;ones(N-2,1)],1) + diag([ ones(N-2,1);2],-1);
        V0=alps;
        tstep =  1;
        t_end = 10;
        tspan = [0:tstep:t_end];
        s0 = V0;
        [T,S] = ode15s(@deRHS,tspan, s0, odeset('maxstep',1));
        u = S(end,:);

        % find the eigenvalues
        Amat=Dscal*A+diag(fp(u));
        onedcrit(j,k) =max(eig(Amat));

        % Now the spherical case
        A =  diag([-2*ones(N-1,1);-2+2*dx/L]) + diag([ ones(N-1,1)],1) + diag([ ones(N-2,1);2],-1);
        [T,S] = ode15s(@deRHSspher,tspan, s0, odeset('maxstep',1));
        u = S(end,:);

        % find the eigenvalues
        Amat=Dscal*A+diag(fp(u./x'));
        threedcrit(j,k)=max(eig(Amat));

    end
end

contour(Alf,Sc,onedcrit',[0 0],'--','linewidth',2)
hold on

contour(Alf,Sc,threedcrit',[0 0],'linewidth',2)
hold off
xlabel('a')
ylabel('\sigma','fontsiz',22)
text(0.05,0.5,'stable','fontsize',20)
text(0.15,2,'unstable','fontsize',20)


end % of main

% ---------------------------------------------------

%%
%the right hand side for ode simulation:
function s_prime=deRHS(t,u)
    global Dscal   A
    s_prime= Dscal*A*u +  f(u);
end
%%
%the right hand side for ode simulation:
function s_prime=deRHSspher(t,u)
    global Dscal   A  x
    s_prime= Dscal*A*u + x.*f(u./x);
end
%%
function out = f(u)
    global gam alps
    w=(u-alps)/gam;
    out = 10*u.*(1-u).*(u-0.5) -w ;
end

%%
function out = fp(u)
    out = 10*(-3.*u.^2 + 3.*u -0.5);
end
