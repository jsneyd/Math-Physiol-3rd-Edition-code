
% Program to solve the time dependent pump-leak model of
% volume control, from Keener and Sneyd, Chapter 2. This is the scaled
% version of the model, so all the parameters (except time) are
% dimensionless.

% Keep in mind that all the concentrations are scaled by Ce, so sometimes
% the values may look a little surprising.

% A couple of functions plot solutions from the full model, obtained by
% solving the 5 ODEs, as described in the book. We check that the correct
% quantity is conserved (which it is, Huzzah!).

% We also plot the steady-state solution obtained by making the assumption
% of intracellular and extracellular electroneutrality. As described in the
% book, this gives an analytic solution for the volume as a function of
% pump rate. When Cm is small, these two solutions agree closely, but
% diverge significantly if Cm is larger (as then the differential equations
% do not preserve electroneutrality, while the algebraic solution does).

function pump_leak

close all; clear all; clc; format longg;

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
par.Ke = 0.06;
par.Nae = 1-par.Ke;     % Choose external solution to be electroneutral
par.z = -1;

par.tau = 1;
par.gamma = 0.11;
par.delta = 0.1;
par.tauw = 1;
par.tauv = 0.001;

options = odeset('RelTol', 1e-11, 'AbsTol', 1e-11);
mu_0       = 1;                                             % cell volume
n_0        = 1/3;                                           % Na in the cell
k_0        = 1;                                             % K in the cell
c_0        = n_0 + k_0 + par.z/mu_0;                        % Determined by electroneutrality. Make sure it's positive!
v_0        = -3;                                            % membrane potential

IC = [mu_0,n_0,k_0,c_0,v_0];                                % The initial condition

singlesolution(IC,par,options)

% If you want to calculate the response to an osmolarity increase,
% uncomment the next line
osmolarity_increase(IC,par,options)

% comparison of numerical and analytic curves
volume_against_pumprate(IC,par,options)

end % of main


%% Plot time-dependent solutions for a single value of P

function singlesolution(IC,par,options)
P = 2;                                                      % Don't make this too big, or the solution breaks.
f_secretion = @(t,x) pump_leak_fns(t,x,par,P);
[t,U] = ode15s(f_secretion, [0,5], IC, options);

plotsolutions(t,U,par)
end

%% Plot time-dependent solutions for an extracellular osmolarity increase

function osmolarity_increase(IC,par,options)
P = 2;                                                      % Don't make this too big, or the solution breaks.
f_secretion = @(t,x) pump_leak_fns(t,x,par,P);
[t1,U1] = ode15s(f_secretion, [0,5], IC, options);

% Now increase the extracellular osmolarity (i.e., increase Ce) and restart the run at the end
% of the previous run. This is tricky, as we are working in scaled
% variables. So the parameters don't change (except for tauw), it's the
% variables that suddenly change.

f = 2 ;                                                     % factor by which we increase Ce
par.tauw = par.tauw/f;
f_secretion = @(t,x) pump_leak_fns(t,x,par,P);

mu_0       = U1(end,1)*f;                                   % cell volume
n_0        = U1(end,2)/f;                                   % Na in the cell
k_0        = U1(end,3)/f;                                   % K in the cell
c_0        = U1(end,4)/f;                                   % Determined by electroneutrality. Make sure it's positive!
v_0        = U1(end,5);                                     % membrane potential

IC = [mu_0,n_0,k_0,c_0,v_0];                                % The new initial condition
[t2,U2] = ode15s(f_secretion, [5,10], IC, options);

t = [t1;t2];
U = [U1;U2];
plotsolutions(t,U,par)
end

%% Plot volume as a function of pump rate

function volume_against_pumprate(IC,par,options)
P = linspace(0.05,8.5,100);  % solution doesn't exist if P is bigger. Be careful with this range.
for i=1:100
f_secretion = @(t,x) pump_leak_fns(t,x,par,P(i));
[t,U] = ode15s(f_secretion, [0,2000], IC, options);
mu(i) = U(end,1);
end

figure(2)
plot(P,mu,'LineWidth',2)
hold on

%--------------------------
% Now add to the graph the curve of volume against pump rate calculated
% analytically by assuming intracellular and extracellular
% electroneutrality, as described in the book.

% If you set up the full model with extracellular electroneutrality, and your initial
% also has intracellular electroneutrality, these two curves are very
% similar, but not identical. However, if you don't have these two conditions,
% the two curves diverge a lot more.

P1 = linspace(0.05,13,1000);
alpha = (par.Nae*exp(-3*P1) + par.Ke*exp(2*P1*par.gamma))/(par.Nae + par.Ke);
a = 4*(1-alpha);                        % quadratic coefficients
b = -4;
c = 1-par.z^2;
mu = (-b + (b^2-4*a*c).^0.5)./(2*a);    % solution of quadratic equation
plot(P1,mu,'r--','LineWidth',2)
xlabel('P')
ylabel('\mu')
ylim([0,10])
xlim([0,10])
set(gca,'FontSize',14)

end


%%  The ODEs

function dx = pump_leak_fns(~,x,par,P)

mu = x(1);
n = x(2);
k = x(3);
c = x(4);
v  = x(5);

VCl = log( 1 / c );                                           % Nernst Potentials
VK = log( par.Ke / k );
VNa = log( par.Nae / n );

pump = P;
ICl = (v + VCl)/par.delta ;
IK = (v - VK)/(par.gamma) - 2*pump ;
INa = v - VNa + 3*pump ;

% Water flux
Q =  n + k + c + 1/mu - 2;

%%% Equations
dx(1) = Q/par.tauw;
dx(2) = (-(INa/par.tau) - dx(1)*n )/mu ;
dx(3) = (-(IK/par.tau) - dx(1)*k)/mu;
dx(4) = (ICl/par.tau - dx(1)*c)/mu;
dx(5) =  ( 1/par.tauv)*( -INa - IK - ICl  );

dx = dx';

end

%% Plot the solutions
function plotsolutions(t,U,par)

mu      = U(:,1);
n 		= U(:,2);
k 		= U(:,3);
c       = U(:,4);
v       = U(:,5);

check = par.tauv*v - par.tau*mu.*(n+k-c);       % Check that this quantity is conserved
electrocheck = n+k-c-1./mu;                     % electroneutrality

figure(1)
subplot(2,4,1)
plot(t,c,'LineWidth',2)
ylabel('Cl')
xlabel('time')
set(gca,'FontSize',14)

subplot(2,4,2)
plot(t,mu,'LineWidth',2)
ylabel('cell volume')
xlabel('time')
set(gca,'FontSize',14)

subplot(2,4,3)
plot(t,v,'LineWidth',2)
ylabel('V')
ylim([-5,0])
xlabel('time')
set(gca,'FontSize',14)

subplot(2,4,4)
plot(t,k,'LineWidth',2)
ylabel('K')
xlabel('time')
set(gca,'FontSize',14)

subplot(2,4,5)
plot(t,n,'LineWidth',2)
ylabel('Na')
xlabel('time')
set(gca,'FontSize',14)

subplot(2,4,6)
plot(t,check,'LineWidth',2)
ylabel('conservation')
xlabel('time')
set(gca,'FontSize',14)

subplot(2,4,7)
plot(t,electrocheck,'LineWidth',2)
ylabel('electroneutrality')
xlabel('time')
set(gca,'FontSize',14)

set(gcf,'position',[900,500,1400,800])
set(gca,'FontSize',14)
end
