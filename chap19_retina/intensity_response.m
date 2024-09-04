
%    -------------------------------------------------------------------
%
%     Intensity-response curves in the model of light adaptation in
%     in cones.
%
%     For Chapter 19, Section 19.2.2 and Exercise 19.3 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

%%
function intensity_response
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);

par.Vstar=35.7;
par.s1=1.59/par.Vstar;  par.s2=1130;  par.vK=-13/par.Vstar;
par.tauy=0.07; par.k1=35.4; par.gam=303; par.delta=5; par.kappa=0.1;
par.eta=52.5; par.tau1=0.012; par.taum=0.02; par.tauz=0.04;

% The IR curves are computed by first holding the cone at a fixed
% background light level (of I0_base) so that the cone reaches steady state
% for this background. Then the background is either increased or decreased
% and the peak response (either the max or the min) is recorded.

% Because we are doing the perturbation using the background light levels,
% stim is always held at zero. This avoids complications due to the kernel
% getting too large, and the light stimulus going negative. Which would be
% bad.

baselist = [0.0001,0.001,0.01,0.1,1];
I0log = linspace(-5,1,100);
vkeep = I0log';

for j=1:5
    I0_base = baselist(j) ;

    I0list = 10.^I0log;

    % First find the steady state
    IC = [0 1 1 1 0];
    tspan = linspace(0,5,200);
    par.I0 = I0_base;
    [t,U] = ode15s(@(t,y)rhs(t,y,par,0), tspan, IC);
    IC = U(end,:);  % initial condition for all subsequent runs

    tspan = linspace(0,2,200);
    for j =1:length(I0list)
        par.I0 = I0list(j);
        [t,U] = ode15s(@(t,y)rhs(t,y,par,0), tspan, IC);
        if I0_base < par.I0
            res(j) = min(U(:,5));
        else
            res(j) = max(U(:,5));
        end
        %plot(t,U(:,5))  % a useful check that the responses look OK
        %hold on
    end
    figure(2)
        plot(I0log,-res*par.Vstar)
        hold on
        plot(log10(I0_base),-par.Vstar*IC(5),'ro') % plot the steady-state value for I0_base
        xlabel('log(I_0)')
        ylabel('-V')
        box off
    % vkeep = [vkeep -(res')*par.Vstar];   % for external plotting

    % Add the Naka-Rushton equation
    nakarushton = 16*I0list./(I0list + 0.015);
    figure(2)
        plot(I0log,nakarushton,'--k','LineWidth',3)
end
%writematrix([vkeep nakarushton'],'IR.dat'); % for external plotting

end % of main




%% The ODEs
function dUdt=rhs(t,U,par,stim)

    p = U(1); x=U(2); y=U(3); z=U(4); v=U(5);
    kernel = (par.eta/par.tau1/6).*((t/par.tau1).^3).*exp(-t/par.tau1);
    s = par.eta*par.I0 + stim*kernel;

    phi = getphi(y,par);

    dUdt(1) = s*(1-p)-par.k1*p;
    dUdt(2) = phi-(par.gam-par.delta)*x*p-par.delta*x;
    dUdt(3) = ((x^3)*exp(-v)-y)/par.tauy;
    dUdt(4) = (((1-par.kappa)/(1+par.kappa))*(x^3)*exp(-v)...
       +(2*par.kappa/(1+par.kappa))*y-z)/par.tauz;
    dUdt(5) = ( (x^3)*exp(-v) - ((1+par.kappa)/3)*z +(par.kappa/2)*y + ...
       ((4+par.kappa)/(6*par.vK))*(v-par.vK) )/par.taum;

    dUdt = dUdt';
end

%%  get phi
function phi = getphi(y,par)
    v = par.vK*(1-y);
    x = (y.*exp(v)).^0.33333;
    I = (exp(-v/par.s1)-1)/par.s2;
    p = par.eta*I./(par.k1+par.eta*I);

    % Comment out the one you don't want
    phi = x.*(par.delta + (par.gam-par.delta)*p);
    %phi = 4 + 84./(1 + (y/0.34).^4);   % Using A(y) instead
end




