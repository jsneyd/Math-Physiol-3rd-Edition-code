%  FHN adjoint solution
function fhn_adjoint_solution
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

 % initial data for integration
 s = [0,0.04];
   
[T,S] = ode15s(@deRHS,tspan,s); 
 
% use the last values to start again, this time stopping when the initial
% values of v and w repeat.

v0=S(end,1);
w0 = S(end,2);

s = [v0,w0];
t=0;
s_prime=deRHS(t,s)
dir = sign(s_prime(1))
stop_cond = odeset('Events',@stopping); 

[Tsv,Sv] = ode15s(@deRHS,tspan,s,stop_cond); 
 
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
  axis([0 period  -0.4 1 ])
  legend('v','w','fontsize',18)
  
% this is the periodic solution
 
% now solve the adjoint equation

s = [0,10]; %initial values
  tstep = -0.001; % integration step size, backwards in time 
t_end =25*period; % length of  interval  a multiple of the period
tspan = [0:tstep:-t_end];   %integrate backwards in time
[Ta,Sa] = ode15s(@deadj_RHS,tspan,s);

% now take the final values and integrate through one more period to find
% the periodic solution

s = [Sa(end,1),Sa(end,2)]; %initial values
tstep = -0.01*period; % integration step size, backwards in time 
t_end =period; % length of  interval  a multiple of the period
tspan = [0:tstep:-period];   %integrate backwards in time
[Ta,Sa] = ode15s(@deadj_RHS,tspan,s);

% figure(2)
% plot(Sa(:,1),Sa(:,2),'linewidth',2)
% xlabel('V','fontsize',20)
% ylabel('W','fontsize',20)
figure(3)
plot(Ta+period,Sa(:,1),Ta+period,Sa(:,2),'linewidth',2)
legend('Va','Wa','fontsize',18)
xlabel('T','fontsize',20)
axis([0 period -30 5])
% check the identity
for j = 1:length(Ta)
   [vt,wt]=  eval_soln(Ta(j))
    FF= deRHS(Ta(j),[vt,wt])
Fvs(j) = FF(1);
Fws(j) = FF(2);
end
id = Fvs'.*Sa(:,1)+Fws'.*Sa(:,2);

figure(4)
plot(Ta+period,id)

% id is supposed to be a constant, but since it is not exactly constant,
% find its average value
idm=mean(id)
figure(4)
plot(Ta+period,id,Ta+period,idm*ones(length(Ta),1),'--','linewidth',2)

%idm is the normalization constant
 
%Define the condition under which to stop the integration 
function [value, isterminal, direction] = stopping(x,y)
global v0 dir
value = [y(1)-v0];  %v0 is the starting value; stop whep the starting value repeats
isterminal = [1];
direction = [dir];
 
function s_prime=deRHS(t,s)  % right hand side for ode system
global alpha eps I
  
 v= s(1);
 w = s(2);

 Fv = (cubic(v)-w)/eps +I;
 Fw = v-w;

s_prime = [Fv; Fw];

function f = cubic(v)
global alpha
f = v.*(v-1).*(alpha-v);
 
function s_prime=deadj_RHS(t,s)  % right hand side for adjoint ode system
global eps 
  
 va= s(1);
 wa = s(2);

 [vt,wt] = eval_soln(t);
 Fva = -cubic_prime(vt)/eps*va-wa; 
 Fwa = va/eps+wa;

s_prime = [Fva; Fwa];

%this evaluates the periodic solution via linear interpolation
function [vt,wt] = eval_soln(t)
global period Tsv Sv
T = mod(t,period);  
tdx = find(Tsv>=T);
ndx = min(tdx);
if (ndx==1)
vt = Sv(1,1);
wt = Sv(1,2);

else
    t0 = Tsv(ndx-1);
    t1 = Tsv(ndx);
    v0 = Sv(ndx-1,1);
    v1 = Sv(ndx,1);
    w0 = Sv(ndx-1,2);
    w1 = Sv(ndx,2);
    vt = v1*(T-t0)/(t1-t0) + v0*(T-t1)/(t0-t1);
    wt = w1*(T-t0)/(t1-t0) + w0*(T-t1)/(t0-t1);
end

function fp = cubic_prime(v)
global alpha
fp = -3*v.^2+2*(alpha+1)*v-alpha;