% code to simulate phase equations
function coupled_oscillators
global N delta Omega

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
N = 32; % number of coupled oscillators
% pick a value of delta
% 
delta = 45;
%delta = 18;
Omega = -10/31;

% set up the differential equation solve

tstep = 1; % integration step size
t_end =2000; % length of  interval
tspan = [0:tstep:t_end];

 % initial data for integration
 s = zeros(1,N);
   
[T,S] = ode15s(@deRHS,tspan,s); 
figure(1)
plot(T,S,'linewidth',2)
xlabel('time','fontsize',20)
ylabel('\phi_i(t)')
formatSpecF = '%6.2f\n';
figure(2)
plot(S(end,:)/T(end),'*','linewidth',2)
xlabel('Oscillator #','fontsize',20)
ylabel('Average Phase Difference')
title(strcat('\delta = ',sprintf(formatSpecF,delta)),'fontsize',18)

function s_prime=deRHS(t,s)  % right hand side for ode system
 global N delta Omega
 K = -2*diag(ones(N,1))+diag(ones(N-1,1),-1) + diag(ones(N-1,1),1);
 
 f = sin(s);
 F = Omega + delta*K*f;


s_prime = F;
