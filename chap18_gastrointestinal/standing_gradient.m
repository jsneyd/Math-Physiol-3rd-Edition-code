
%    -------------------------------------------------------------------
%
%     Solve the standing gradient water transport equations via
%     shooting
%
%     For Chapter 18, Section 18.3 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function standing_gradient
close all
clear all
clc
global D P r c0 alp L N0
set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0, ...
'defaultpatchlinewidth', 0.7);

%parameters
D=1000;
r=0.05;
c0 = 0.3;
P=0.2;
L=100;
alp=0.1;

N0list = [0.1,0.3];
for j = 1:length(N0list)
    N0 = N0list(j);
    % guess a lower bound for the unknown initial value
    a = 0.1;
    % guess an upper bound for the unknown initial value
    b = 2;

    %call the function bisect(a,b,@fval)
    tic
    root = bisect(a,b,@transport_eqs);
    toc
    % now plot the solution
    [X,S] = desolve(root);
    Nb = alp*L;
    % Using plotyy is deprecated Matlab, but necessary for Octave, which
    % doesn't yet work with yyaxis. If you are running this in Matlab, not
    % Octave, you really should use plotyy instead. It's much better.
    figure(2)
        [hAx,hLine1,hLine2] = plotyy(X,S(:,2),X,S(:,1));
        ylabel(hAx(2),'v (mm/s)')
        ylabel(hAx(1),'c(x) (mM)','fontsize',20)
        xlabel('x (mm)','fontsize',20)
        ylim(hAx(1),[0.2,1.3])
        ylim(hAx(2),[0,300])
        box off
        hold on
end
hold off

%Now solve the BVP for variable L
Llist = [10:0.501:100];
for j = 1:length(Llist)
    L = Llist(j);
    N0 = 0.1;

    root = bisect(a,b,@transport_eqs);
    [X,S] = desolve(root);
    ce(j) = S(end,2)-D*S(end,3)/S(end,1);
end
figure(3)
    plot(Llist,ce)
    xlabel('L (mm)')
    ylabel('c_e (mM)')
    axis([10 100 0 3])
box off
end % of main


%%
% the function transport_eqs(c) solves the initial value problem with
% initial value c and calculates the error for the boundary value condition
% at x=L
function fval = transport_eqs(c)
    global D P r c0 alp L N0

    [X,S] = desolve(c);
    fval = S(end,2)-c0;
end

%%
function [X,S] = desolve(c)
    global D P r c0 alp L N0
    xstep = L/100; % integration step size
    x_end = L; % length of  interval
    xspan = [0:xstep:x_end];

    % initial data for integration
    y0 = [0;c;0];

    stop_cond = odeset('Events',@stopping);   % The stopping conditions for the integration
    %this is needed since for some initial conditions, the integrator cannot
    %integrate all the way to x=L
    [X,S] = ode23(@deRHS,xspan,y0,stop_cond);
    % to save time, don't plot these solutions
    % figure(1)
    % plot(X,S(:,2))
    % hold on
end

%% specify the right hand side of the differential equation system
function s_prime=deRHS(x,s)  % right hand side for ode system
    global D P r c0 alp L N0
    Nx = N0*(x<=alp*L);
    v=s(1);
    c = s(2);
    cp = s(3);

    vp = 2*P*(c-c0)/r;
    Fc = cp;
    Fcp = (vp*c+v*cp-2*Nx/r)/D;
    s_prime = [vp,Fc,Fcp]';
end

%% Define the condition under which to stop the integration
function [value, isterminal, direction] = stopping(x,y)
    value = [y(2);y(2)-5];
    isterminal = [1;1];
    direction = [0;0];
end

%% this is the bisection algorithm
function root = bisect(a,b,feval)
    ul = a;
    fl = feval(ul);
    uu = b;
    fu = feval(uu);

    % we make the assumption, without checking, that
    % fu*fl<0

    % if not, the algorithm fails to find a root.
    N = 45;  % number of iterates
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



