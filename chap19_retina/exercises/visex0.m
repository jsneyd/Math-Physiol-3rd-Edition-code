%%
clear all
close all
clc

par.I0=0;
stims = [0.001 0.01 0.1 1];
par.s1=1.59;  par.s2=1130;  par.E=-13;  par.Vstar=35.7;
par.tauy=0.07; par.k1=35.4; par.gam=303; par.delta=5; par.kappa=0.1;
par.eta=52.5; par.tau1=0.012; par.taum=0.016; par.tauz=0.04;
par.V0=-par.s1*log(1+par.s2*par.I0);

% First find the steady state
tspan = linspace(0,2,2000);
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
figure(1)
plot(t,U(:,5))

%% The ODEs

function dUdt=rhs(t,U,par,stim)

p = U(1); x=U(2); y=U(3); z=U(4); V=U(5);

kernel = (par.eta/par.tau1/6).*((t/par.tau1).^3).*exp(-t/par.tau1);
s = par.eta*par.I0 + stim*kernel;

phi=(y*exp(par.E*(1-y)/par.Vstar))^(1/3)* ...
        (par.delta+(par.gam-par.delta)*par.eta* ...
        (exp(-par.E*(1-y)/par.s1)-1)/par.s2/ ...
        (par.k1+par.eta*(exp(-par.E*(1-y)/par.s1)-1)/par.s2));

dUdt(1) = s*(1-p)-par.k1*p;
dUdt(2) = phi-(par.gam-par.delta)*x*p-par.delta*x;
dUdt(3) = ((x^3)*exp(-V/par.Vstar)-y)/par.tauy;
dUdt(4) = (((1-par.kappa)/(1+par.kappa))*(x^3)*exp(-V/par.Vstar)...
   +(2*par.kappa/(1+par.kappa))*y-z)/par.tauz;
dUdt(5) = ((-6*par.E/(4+par.kappa))*(x^3)*exp(-V/par.Vstar)+ ...
   (2*(1+par.kappa)/(4+par.kappa))*par.E*z-(V-par.E)-...
   (3*par.E*par.kappa/(4+par.kappa))*y)/par.taum;

dUdt = dUdt';

end


