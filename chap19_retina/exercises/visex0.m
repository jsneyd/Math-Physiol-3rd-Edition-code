%%
clear all
close all
clc

par.Vstar=35.7;
par.s1=1.59/par.Vstar;  par.s2=1130;  par.vK=-13/par.Vstar;
par.tauy=0.07; par.k1=35.4; par.gam=303; par.delta=5; par.kappa=0.1;
par.eta=52.5; par.tau1=0.012; par.taum=0.02; par.tauz=0.04;

%plotphi(par)
%plot_impulse_sequence(par)
plot_single_impulse_responses(par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the function phi -------------------------------------------------
function plotphi(par)

y = linspace(0.2,1,500);
phi = getphi(y,par);
A = 4 + 84./(1 + (y/0.34).^4);
figure(1)
plot(y,phi,y,A,'r--','LineWidth',2)
legend('\phi(y)', 'A(y)')
set(gca,'FontSize',14)
xlabel('y')

end

%% Plot a sequence of impulse responses -------------------------------------------------
function plot_impulse_sequence(par)
par.I0=0.01;

% First find the steady state
tspan = linspace(0,2,200);
IC = [0 1 1 1 0];
[t1,U1] = ode15s(@(t,y)rhs(t,y,par,0), tspan, IC);
IC = U1(end,:);
[t2,U2] = ode15s(@(t,y)rhs(t,y,par,0.0001), tspan, IC);
IC = U2(end,:);
[t3,U3] = ode15s(@(t,y)rhs(t,y,par,0.001), tspan, IC);
IC = U3(end,:);
[t4,U4] = ode15s(@(t,y)rhs(t,y,par,0.01), tspan, IC);
IC = U4(end,:);
[t5,U5] = ode15s(@(t,y)rhs(t,y,par,0.1), tspan, IC);

t = [t1; t2+2; t3+4; t4+6; t5+8];
U = [U1;U2;U3;U4;U5];
figure(2)
plot(t,U(:,5))

end

%% Plot some single impulse responses -------------------------------------------------
function plot_single_impulse_responses(par)

par.I0=0.0005;
% First find the steady state
tspan = linspace(0,2,500);
IC = [0 1 1 1 0];
[t1,U1] = ode15s(@(t,y)rhs(t,y,par,0), tspan, IC);
% Then add the impulse response
IC = U1(end,:);
[t2,U2] = ode15s(@(t,y)rhs(t,y,par,par.I0), tspan, IC);

figure(3)
plot(t2,U2(:,5)-IC(5),'LineWidth',2)
%plot(t2,(U2(:,5)-IC(5))/max(abs(U2(:,5)-IC(5))),'LineWidth',2)
hold on


par.I0=0.005;
% First find the steady state
tspan = linspace(0,2,500);
IC = [0 1 1 1 0];
[t1,U1] = ode15s(@(t,y)rhs(t,y,par,0), tspan, IC);
% Then add the impulse response
IC = U1(end,:);
[t2,U2] = ode15s(@(t,y)rhs(t,y,par,par.I0), tspan, IC);
plot(t2,U2(:,5)-IC(5),'LineWidth',2)
%plot(t2,(U2(:,5)-IC(5))/max(abs(U2(:,5)-IC(5))),'LineWidth',2)

par.I0=0.05;
% First find the steady state
tspan = linspace(0,2,500);
IC = [0 1 1 1 0];
[t1,U1] = ode15s(@(t,y)rhs(t,y,par,0), tspan, IC);
% Then add the impulse response
IC = U1(end,:);
[t2,U2] = ode15s(@(t,y)rhs(t,y,par,par.I0), tspan, IC);
plot(t2,U2(:,5)-IC(5),'LineWidth',2)
%plot(t2,(U2(:,5)-IC(5))/max(abs(U2(:,5)-IC(5))),'LineWidth',2)

par.I0=0.5;
% First find the steady state
tspan = linspace(0,2,500);
IC = [0 1 1 1 0];
[t1,U1] = ode15s(@(t,y)rhs(t,y,par,0), tspan, IC);
% Then add the impulse response
IC = U1(end,:);
[t2,U2] = ode15s(@(t,y)rhs(t,y,par,par.I0), tspan, IC);
plot(t2,U2(:,5)-IC(5),'LineWidth',2)
%plot(t2,(U2(:,5)-IC(5))/max(abs(U2(:,5)-IC(5))),'LineWidth',2)

xlim([0 0.5])
xlabel('seconds')
ylabel('v')
set(gca,'FontSize',14)
legend('I_0=0.0005','I_0=0.005','I_0=0.05','I_0=0.5','Location','SouthEast')
box off

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
%phi = x.*(par.delta + (par.gam-par.delta)*p);
phi = 4 + 84./(1 + (y/0.34).^4);   % Using A(y) instead, just for fun
end


