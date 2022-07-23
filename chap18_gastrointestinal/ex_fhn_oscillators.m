% coupled FHN oscillators
function fhn_oscillators
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

global N   alpha G delj K

N = 200; % number of coupled oscillators
Nj = [1:N]';
delj = 0.02+0.06*Nj/N;
% parameters
 alpha = 0.4;
 % coupling coefficient

Glist = [0;50];
tlengthlist = [50,300];
  %K is the coupling matrix
 K = - diag([1; 2*ones(N-2,1);1])+diag(ones(N-1,1),-1) + diag(ones(N-1,1),1);
% sketch a phase plane:
v =[-0.2:.01:1.2];
f = cubic(v);
figure(1)
plot(v,f,[alpha,alpha],[-0.2,0.3],'linewidth',2)
hold on

for gj = 1:length(Glist)
    G = Glist(gj);
% set up the differential equation solve

tstep = 0.01; % integration step size
t_end =tlengthlist(gj); % length of  interval
tspan = [0:tstep:t_end];

 % initial data for integration
 s = zeros(1,2*N);
   
[T,S] = ode15s(@deRHS,tspan,s); 
figure(1)
plot(S(:,N),S(:,2*N),'--',S(:,1),S(:,N+1),'--')
xlabel('v','fontsize',20)
ylabel('w','fontsize',20)
 
% find the oscillation periods
for j = 1:N
Q = (S(2:end,j)-alpha).*(S(1:end-1,j)-alpha);
%Check when v crosses alpha
Tj = find(Q<=0);
 
% calculate the frequency
nj = length(Tj);

np =fix(nj/2);
timetot = T(Tj(np)) - T(Tj(2)); 
per(j) = timetot/(np-1);
end

freq = 1./per;
figure(4)
 
if(gj==1)
plot([1:N],freq,'linewidth',2)
else
plot([1:N],freq,'*') 
end
hold on
xlabel('oscillator #','fontsize',20)
ylabel('frequency','fontsize',20)
 text(100,1.45,'coupled','fontsize',18)
 text(50,1.3,'uncoupled','fontsize',18)

 
end
function s_prime=deRHS(t,s)  % right hand side for ode system
 global  N alpha G delj K
  
 v= s(1:N);
 w = s(N+1:2*N);

 
 Fv = (cubic(v)-w)./delj +G*K*v;
 Fw = v-alpha;

s_prime = [Fv; Fw];

function f = cubic(v)
global alpha
f = v.*(v-1).*(alpha-v);
