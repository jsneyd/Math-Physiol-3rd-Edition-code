function fhn_periodic_solution_w_Floquet
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

global alpha eps I v0 Sv Tsv period dir

% parameters
 alpha = 0.05;
 eps = 0.01;
 I = 5;

% first make a phase portrait
vmin = -0.4;
vmx = 1.1;

v =[vmin:.01:vmx];
f = cubic(v);
ff = f+eps*I;
figure(1)
plot(v,ff,v,v,'linewidth',2)

hold on
% set up the differential equation solve

tstep = 0.01; % integration step size
t_end =100; % length of  interval; run for a long time
tspan = [0:tstep:t_end];

 % initial data for integration  run for a long time to converge to
 % periodic solution
 s = [0,0.04];
   
[T,S] = ode15s(@deRHS,tspan,s); 
 
% use the last values to start again,  stopping when the initial
% values of v and w repeat.

v0=S(end,1);
w0 = S(end,2);

s0 = [v0,w0]
t=0;
s_prime=deRHS(t,s0)
dir = sign(s_prime(1))
stop_cond = odeset('Events',@stopping); 
s = [v0,w0,1,0,0,1];  % include initial data for fundamental solution
[Tsv,Sv] = ode15s(@deRHS_lin,tspan,s,stop_cond); 
 
 period = Tsv(end);
  
 figure(1)
 plot(Sv(:,1),Sv(:,2),'linewidth',2)
 xlabel('v','fontsize',20)
 ylabel('w','fontsize',20)
 title('Periodic solution - phase portrait','fontsize',20)
 axis([vmin vmx 0 max(f)])
 hold off

 figure(2)
 plot(Tsv,Sv(:,1),Tsv,Sv(:,2),'linewidth',2)
axis([0 period  -0.4 1])
xlabel('t','fontsize',20)
legend('v','w','fontsize',20)
title('Periodic solution','fontsize',20)

M=[Sv(end,3),Sv(end,4);Sv(end,5),Sv(end,6)];
period = Tsv(end)
E=eig(M);
E(1)
E(2)
 function s_prime=deRHS(t,s)  % right hand side for ode system
global eps I
  
 v= s(1);
 w = s(2);

 Fv = (cubic(v)-w)/eps +I;
 Fw = v-w;
s_prime = [Fv;Fw];
function s_prime=deRHS_lin(t,s)  % right hand side for ode system with linearization
global eps I
 
 v= s(1);
 w = s(2);
 V1 = s(3);
 W1 = s(4);
 V2 = s(5);
 W2 = s(6);

 Fv = (cubic(v)-w)/eps +I;
 Fw = v-w;

 FV1 = (cubic_prime(v)*V1-W1)/eps;
 FW1 = V1-W1;
 FV2 = (cubic_prime(v)*V2-W2)/eps;
 FW2 = V2-W2;
 
s_prime = [Fv;Fw;FV1;FW1;FV2;FW2];

function f = cubic(v)
global alpha
f = v.*(v-1).*(alpha-v);

function fp = cubic_prime(v)
global alpha
fp = -3*v.^2+2*(alpha+1)*v-alpha;
 
 %Define the condition under which to stop the integration 
function [value, isterminal, direction] = stopping(x,y)
global v0 dir
value = [y(1)-v0];  %v0 is the starting value; stop whep the starting value repeats
isterminal = [1];
direction = [dir];
 