
function phase_splitting
global eps
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

%parameters

eps = 0.01;

phi = [0:.01:4]*pi;

rhs = -2-eps-2*sin(phi);
figure(1)
plot(phi/pi,rhs,'linewidth',2)
xlabel('\phi/\pi','fontsize',20)
ylabel('d\phi/dt','fontsize',20)
tstep = 0.01; % integration step size
t_end =100; % length of  interval
tspan = [0:tstep:t_end];

 % initial data for integration
 s = 0;
   
[T,S] = ode15s(@deRHS,tspan,s); 

figure(2)
plot(T,S/pi,'linewidth',2)
xlabel('time','fontsize',20)
ylabel('\phi/\pi','fontsize',20)

function s_prime=deRHS(t,s)  % right hand side for ode system
 global  eps
F = -2-eps-2*sin(s)
s_prime = F;