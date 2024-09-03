 
%    -------------------------------------------------------------------
%
%     Model of light adaptation in cones.
%
%     For Chapter 19, Section 19.2.2 and Exercise 19.3 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

%%
function adaptation
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

plot_impulse_sequence(par)

end % of main


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot a sequence of impulse responses -------------------------------------------------
function plot_impulse_sequence(par)
 formatSpecF = '%6.2f\n';
    par.I0=0;
    I0list = [0,0.1,1];
    for k = 1:length(I0list)
        par.I0 = I0list(k)
    stimlist = [0,0.0001,0.001,0.01,0.1];
    U=[];
    t=[];
    IC = [0 1 1 1 0];
    for j =1:length(stimlist)
        stim = stimlist(j);
        tspan = linspace(0,2,200);       
        [t1,U1] = ode15s(@(t,y)rhs(t,y,par,stim), tspan, IC);
        U = [U;U1];
        t = [t;t1+2*(j-1)];   
        IC = U1(end,:);
    end  
    figure(k)
        plot(t,U(:,5)*par.Vstar)
        xlabel('t (s)')
        ylabel('V-V_d')
         title(strcat('I_0 = ',sprintf(formatSpecF,par.I0)),'fontsize',18)
        box off
    end
    %writematrix([t U(:,5)*par.Vstar],'adaptation.dat')
end


%% The ODEs -------------------------------------------------
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

%%  get phi -------------------------------------------------
function phi = getphi(y,par)
    v = par.vK*(1-y);
    x = (y.*exp(v)).^0.33333;
    I = (exp(-v/par.s1)-1)/par.s2;
    p = par.eta*I./(par.k1+par.eta*I);
    
    % Comment out the one you don't want
    phi = x.*(par.delta + (par.gam-par.delta)*p);
    %phi = 4 + 84./(1 + (y/0.34).^4);   % Using A(y) instead 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



