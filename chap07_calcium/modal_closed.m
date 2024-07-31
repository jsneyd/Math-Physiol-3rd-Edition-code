clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 

global ct p tau_max  Ktau  Kc Kh kf Vserca kserca Kbar Kp gamma
global direction   % this parameter is a cheat, so that we can solve the ode backwards

% Parameter values
ct = 2;
p = 0.35;
tau_max=1000;

Ktau=0.1;
tauP=1;
Kc=0.2;
Kh=0.08;
kf=10;
 
Vserca=0.9;
kserca=0.2;
Kbar=0.00001957;
Kp=0.2;

gamma=5.5;

% solve the ode
direction = 1;  % solve forwards in time to get the stable limit cycle.
init = [0.05,0.004];
tspan = linspace(0,150,20000);
[T,Y] = ode15s(@(t,x)oderhs(t,x),tspan,init);
cc1 = Y(end/2:end,1);
hh1 = Y(end/2:end,2);

% plot solution
figure(1)
plot(T(end/2:end),cc1)
hold on

% Now solve backwards in time to find the unstable limit cycle
direction = -1;
init = [0.34,0.0065];
tspan = linspace(0,50,20000);
[T,Y] = ode15s(@(t,x)oderhs(t,x),tspan,init);
cc2 = Y(end/2:end,1);
hh2 = Y(end/2:end,2);


%%
% calculate and plot the nullclines

c = linspace(0.01,1,1000);
ce=gamma*(ct-c);

phi_c=c.^4./(c.^4+Kc^4);
phi_p=p^2/(Kp^2+p^2);
phi_p_down=Kp^2/(Kp^2+p^2);
h_inf=Kh^4./(Kh^4+c.^4);
tauh = tau_max*Ktau^4./(Ktau^4+c.^4);
alpha = phi_p_down.*(1-phi_c.*h_inf);

serca=Vserca*(c.*c-Kbar*ce.*ce)./(c.*c+kserca*kserca);
dum = serca./(kf*(ce-c));
betadum = 0.4*alpha.*dum./(1-dum*(1+0.4));
c_nullcline = betadum./(phi_p*phi_c);
h_nullcline = h_inf;

figure(2)
plot(c,c_nullcline,c,h_nullcline)
xlim([0,1]); xlabel('c');
ylim([0,0.02]); ylabel('h');

% add the solutions to the phase plane
hold on
plot(cc1,hh1,cc2,hh2,'--')
legend('dc/dt=0','dh/dt=0','stable limit cycle','unstable limit cycle')

% save for proper plotting
% igorsave1 = [c' c_nullcline' h_nullcline'];
% igorsave2 = [cc1 hh1 cc2 hh2];
% writematrix(igorsave1,'modal_out1.dat')
% writematrix(igorsave2,'modal_out2.dat')


%% Define the differential equations here

function out=oderhs(t,x)
c = x(1);
h = x(2);

global ct p tau_max  Ktau  Kc Kh kf Vserca kserca Kbar Kp gamma
global direction

ce=gamma*(ct-c);

phi_c=c^4/(c^4+Kc^4);
phi_p=p^2/(Kp^2+p^2);
phi_p_down=Kp^2/(Kp^2+p^2);
h_inf=Kh^4/(Kh^4+c^4);
tauh = tau_max*Ktau^4/(Ktau^4+c^4);

beta = phi_p*phi_c*h;
alpha = phi_p_down*(1-phi_c*h_inf);
Po=beta/(beta + 0.4*(beta+alpha));

serca=Vserca*(c*c-Kbar*ce*ce)/(c*c+kserca*kserca);

dcdt  = kf*Po*(ce-c)-serca;
dhdt  = (h_inf-h)/tauh;
out = direction*[dcdt,dhdt]';

end

