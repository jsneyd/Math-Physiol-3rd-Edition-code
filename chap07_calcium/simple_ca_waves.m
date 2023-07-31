% code to simulate calcium waves with diffusing IP3
% use method of lines


function Diff_sim
 global sc   
 
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
    
 %parameters for calcium dynamics

p.km = 1;
p.ks = 20;
p.a0 = 0.1;
p.a1 = 0.1;
p.ph1 = 2;
p.ph2 = 1;
p.a = 0.05;
p.kf = 20;
p.del = 1;
p.gam = 5;
p.Dc = 25;  %Ca diffusion coefficient
p.p =3; %the bifurcation parameter

% make a phase portrait
c = [0.01:.01:1.5];
Po = p.p*(c.^2./(c.^2+p.ph1^2)).*(p.ph2./(c+p.ph2)); %open probability

 
Jserca = p.ks*c;
Jpm = p.km*c;
Jin = p.a0+p.a1*p.p;

ce1=(Jserca-p.del*(Jin-Jpm))./(p.kf*Po+p.a)+c;
ce2= Jserca./(p.kf*Po+p.a)+c;

% integrate the ode
 init = [0.5 , 3]; %initial data for the ode
tstep = 0.01;
t_end = 50;
 
%specify the output points
tspan = [0:tstep:t_end];
p.N=1;
 
[T,S] = ode15s(@(t,x)coscrhs(t,x,p),tspan,init);
C = S(:,1);
Ce = S(:,2);
 figure(1)
plot(c,ce1,'--',c,ce2,'--',C,Ce)
legend('dc/dt=0','dc_e/dt=0')
xlabel('c')
ylabel('c_e')
axis([0 1.5 0 15])

figure(2)
 plot(T,C,T,Ce,'linewidth',2)
 
xlabel('Time','fontsize',20)
legend('Ca^{++}','Ca^{++}_e')
  
% now integrate the pde
  p.N=100;  % number of spatial grid points
  p.L=5;
  p.h=p.L/p.N;
   
  p.dv =  0.0; % diffusion coefficient for Ca

 p.sc = [1;2*ones(p.N-2,1);1];
  X = p.h*(1:p.N)';
  %set initial data
    
 U = 4*exp(-3*X.^2/p.L);
 V =0.05*ones(p.N,1); 
 y0=1.5*ones(p.N,1);
ce0=10*ones(p.N,1);
  %plot the solution at each time step

 figure(3)
  subplot(2,1,1)
  plot(X,V,'--','linewidth',2)
  xlabel('x','fontsize',20)
   ylabel('v(x,t)','fontsize',20)
  hold on
  subplot(2,1,2)
  plot(X, U,'--','linewidth',2)
   xlabel('x','fontsize',20)
    ylabel('u(x,t)','fontsize',20)
  %axis([0 L  0 0.5])
 % axis([0 L 0 1])
 hold on 
 
tstep = 1; % time between plots
t_end = 40; %total time to run simulation
 
%specify the output points
tspan = [0:tstep:t_end];

% integrate the system of ode's:
[T,S] = ode23( @(t,x)pdeRHS(t,x,p),tspan,[U;V;y0;ce0], odeset('maxstep',1));  

% plot the solution at each time step
 
for j = 1:length(T)
   figure(3)
  subplot(2,1,1)
  plot(X,S(j,p.N+1:2*p.N),'linewidth',2)
  xlabel('x','fontsize',20)
   ylabel('Ca^{++}(x,t)','fontsize',20)
  hold on
  subplot(2,1,2)
  plot(X, S(j,1:p.N),'linewidth',2)
   xlabel('x','fontsize',20)
    ylabel('IP_3(x,t)','fontsize',20)
    
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

  mesh(S(:,p.N+1:2*p.N))

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
 FP = scu*(-p.sc.*P+[0;P(1:end-1)]+[P(2:end);0])%+Fp;
FC = scv*(-p.sc.*C+[0;C(1:end-1)]+[C(2:end);0]) +Fc ;
 

s_prime = [FP;FC;Fy;Fce];
 

%the right hand side for ode simulation:
 
  
function out=coscrhs(t,s,p)
 % evaluate the ode dynamics
% there are two variables
 
c = s(1: p.N); % calcium

ce = s(p.N+1:2*p.N); %ER calcium

Po = p.p*(c.^2./(c.^2+p.ph1^2)).*(p.ph2./(c+p.ph2)); %open probability

Jipr = (p.kf*Po+p.a).*(ce-c);
Jserca = p.ks*c;
Jpm = p.km*c;
Jin = p.a0+p.a1*p.p;

Fc = Jipr-Jserca+p.del*(Jin-Jpm);
 
 
Fce = p.gam*(Jserca-Jipr);
out = [Fc;Fce];
 