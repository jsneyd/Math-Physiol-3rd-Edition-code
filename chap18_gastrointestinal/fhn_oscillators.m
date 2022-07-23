% coupled FHN oscillators
function fhn_oscillators
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

global N delta alpha G delj K

N = 200; % number of coupled oscillators
nj = [1:N]';
delj = 0.02+0.06*nj/N;
% parameters

delta = 0.05;
alpha = 0.4;

  %K is the coupling matrix

 K = - diag([1; 2*ones(N-2,1);1])+diag(ones(N-1,1),-1) + diag(ones(N-1,1),1);
% sketch a phase plane:
v =[-0.2:.01:1.2];
f = cubic(v);
figure(1)
plot(v,f,[alpha,alpha],[-0.1,0.1],'linewidth',2)
hold on

G = 50; % coupling coefficient
 
% set up the differential equation solve

tstep = 0.01; % integration step size
t_end =400; % length of  interval
tspan = [0:tstep:t_end];

 % initial data for integration
 s = zeros(1,2*N);
   
[T,S] = ode15s(@deRHS,tspan,s); 
figure(1)
plot(S(:,N),S(:,2*N),'--',S(:,1),S(:,N+1),'--')
hold off
 
% find the oscillation periods
for j = 1:N
Q = (S(2:end,j)-alpha).*(S(1:end-1,j)-alpha);
%Check when v crosses alpha
Tj = find(Q<=0);
%figure(3)
%plot(T(Tj),j*ones(length(Tj),1),'*')
%hold on

% calculate the frequency
nj = length(Tj);

np =fix(nj/2);
timetot = T(Tj(np)) - T(Tj(2)); 
per = timetot/(np-1);
freq = 1/per;
figure(4)
plot(j,per,'*')
hold on
xlabel('oscillator #','fontsize',20)
ylabel('frequency','fontsize',20)



end

function s_prime=deRHS(t,s)  % right hand side for ode system
 global N delta alpha G delj K
  
 v= s(1:N);
 w = s(N+1:2*N);

 
 Fv = (cubic(v)-w)./delj +G*K*v;
 Fw = v-alpha;

s_prime = [Fv; Fw];
function f = cubic(v)
global alpha
f = v.*(v-1).*(alpha-v);
