% code to simulate calcium waves with diffusing IP3
% use method of lines


function Diff_sim
 global sc   
 
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
    
 %parameters for calcium dynamics

p.Vplc = 0.26;
p.ct = 2;
p.k1 = 400;
p.k2 = 0.2;
p.k3 = 400;
p.k4 = 0.2;
p.k5 = 20;
p.km1 = 52;
p.km2 = 0.21;
p.km3 = 377.2;
p.km4 = 0.0289;
p.km5 = 1.64;
p.K1=p.km1/p.k1;
p.K5=p.km5/p.k5;
p.K2 = p.km2/p.k2;
p.K3=p.km3/p.k3;
p.K4=p.km4/p.k4;

p.gm = 5.5;
p.n = 2;
p.Kbar = 0;
p.kf=1.11;
p.Ki=0.4;
p.kf2=0.0203;
 
p.k3K=0;
p.k5P=0.66;
p.Kplc=0.2;
p.Kplc = 0;  % gives the no feedback case
p.Vserca = 0.9; 
p.Kserca = 0.1;
 
p.tauy=12.5;
 
% integrate the ode
 init = [0.35 ,0.2,1.,9]; %initial data for the ode
tstep = 1;
t_end = 200;
 
%specify the output points
tspan = [0:tstep:t_end];
p.N=1;
 
[T,S] = ode15s(@(t,x)coscrhs(t,x,p),tspan,init);
P = S(:,1);
C = S(:,2);
y = S(:,3);
Ce = S(:,4);
 
figure(1)
 plot(T,P,T,C,T,Ce,T,y,'linewidth',2)
 
xlabel('Time','fontsize',20)
legend('IP_3','Ca^{++}')
  
% now integrate the pde
  p.N=100;  % number of spatial grid points
  p.L=5;
  p.h=p.L/p.N;
  p.du = 0.05; % diffusion coefficient for IP3
  p.dv =  0.0; % diffusion coefficient for Ca

 p.sc = [1;2*ones(p.N-2,1);1];
  X = p.h*(1:p.N)';
  %set initial data
    ipmx = 4;
    U = P(end)*ones(p.N,1)+ipmx*exp(-3*X.^2);  % initial ip3 concentration
 V =C(end)*ones(p.N,1); % ca concentration
 y0=y(end)*ones(p.N,1);
ce0=Ce(end)*ones(p.N,1);
init = [U;V;y0;ce0];
size(init)
  %plot the solution at each time step

 figure(3)
  subplot(2,1,1)
  plot(X,V,'--','linewidth',2)
  xlabel('x','fontsize',20)
   ylabel('v(x,t)','fontsize',20)
   
  subplot(2,1,2)
  plot(X, U,'--','linewidth',2)
   xlabel('x','fontsize',20)
    ylabel('u(x,t)','fontsize',20)
  %axis([0 L  0 0.5])
 % axis([0 L 0 1])
 
 
tstep = 0.5; % time between plots
t_end = 150; %total time to run simulation
 
%specify the output points
tspan = [0:tstep:t_end];

% integrate the system of ode's:
[T,S] = ode23( @(t,x)pdeRHS(t,x,p),tspan,init, odeset('maxstep',1));  

% plot the solution at each time step
 
for j = 1:length(T)
   figure(3)
  subplot(2,1,1)
  plot(X,S(j,p.N+1:2*p.N),'linewidth',2)
  xlabel('x','fontsize',20)
   ylabel('Ca^{++}(x,t)','fontsize',20)
   axis([0 max(X) 0 0.5])
  subplot(2,1,2)
  plot(X, S(j,1:p.N),'linewidth',2)
   xlabel('x','fontsize',20)
    ylabel('IP_3(x,t)','fontsize',20)
     axis([0 max(X) 0 ipmx])
end
 
 figure(5) %plot the final solution
 plot(X, S(end,1:p.N),X,S(end,p.N+1:2*p.N),'linewidth',2)
 legend('IP_3','Ca^{++}','fontsize',20)
 xlabel('x','fontsize',20)

 figure(6)
 %plot a time course
 plot(T,S(:,3*p.N/2))
  xlabel('time')
ylabel('Ca^{++}(L/2)')

figure(7)
  contour(X,T,S(:,p.N+1:2*p.N),'linewidth',2)
  xlabel('x')
  ylabel('t')
  zlabel('Ca^{++}')

  figure(8)
  mesh(X,T,S(:,1:p.N),'linewidth',2)
  xlabel('x')
  ylabel('t')
  zlabel('IP_3')


%the right hand side for pde (MoL) simulation:
function s_prime=pdeRHS(t,s,p)
 
%There are four variables:  P, C, y, ce 
% only P and C are diffusing
p.h = p.L/p.N;
 
scu = p.du/p.h^2;
scv = p.dv/p.h^2;
 
P = s(1:p.N);
C = s(p.N+1:2*p.N);
y = s(2*p.N+1:3*p.N);
ce = s(3*p.N+1:4*p.N);

% this is method of lines
% evaluate the ode part
%[Fp;Fc;Fy;Fce]=
out=coscrhs(t,s,p);
 
Fp=out(1:p.N);
Fc=out(p.N+1:2*p.N);
Fy=out(2*p.N+1:3*p.N);
Fce=out(3*p.N+1:4*p.N);
FP = scu*(-p.sc.*P+[0;P(1:end-1)]+[P(2:end);0]);%+Fp; % no IP3 feedback
FC = scv*(-p.sc.*C+[0;C(1:end-1)]+[C(2:end);0]) +Fc ;
 
s_prime = [FP;FC;Fy;Fce];
 
%the right hand side for ode simulation:
   
function out=coscrhs(t,s,p)
 % evaluate the ode dynamics
% this is  a simple closed cell model
IP3 = s(1:p.N); %IP3
c = s(p.N+1:2*p.N); % calcium
y = s(2*p.N+1:3*p.N);  % inactivation
ce = s(3*p.N+1:4*p.N); %ER calcium

Po = (IP3.*c.*(1-y)./((IP3+p.K1).*(c+p.K5))).^3; %open probability

Jipr = (p.kf*Po+p.kf2).*(ce-c);
Jserca = p.Vserca.*(c.^2-p.Kbar*ce.^2)./(p.Kserca^2+c.^2);

Fp =  p.Vplc*c.^p.n./(p.Kplc.^2+c.^p.n)-(p.k3K+p.k5P).*IP3; %IP3 feedback dynamcs
Fp=zeros(p.N,1);  % IP3 is not reacting
Fc = Jipr-Jserca;
Fy =  (-1+(1-y).*(p.Ki+c)/p.Ki)/p.tauy;
 
Fce = -p.gm*Fc;
out = [Fp;Fc;Fy;Fce];
 