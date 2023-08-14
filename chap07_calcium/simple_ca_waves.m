% code to simulate calcium waves (with diffusing IP3) for a simple 2 variable calcium
% model
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
parlist=[2.5,0.5];%the bifurcation parameter
c0list = [0.5,0.28];
ce0list=[5.2,30];
axlist=[15,35];
for j =1:length(parlist)
    p.p=parlist(j);%the bifurcation parameter

c0=c0list(j);
ce0=ce0list(j);

% make a phase portrait for each parameter on the list
c = [0.01:.01:1.5];
Po = p.p*(c.^2./(c.^2+p.ph1^2)).*(p.ph2./(c+p.ph2)); %open probability

 
Jserca = p.ks*c;
Jpm = p.km*c;
Jin = p.a0+p.a1*p.p;

ce1=(Jserca-p.del*(Jin-Jpm))./(p.kf*Po+p.a)+c;
ce2= Jserca./(p.kf*Po+p.a)+c;
 
% integrate the ode
 init = [c0,ce0]; %initial data for the ode
 ct = c0+ce0/p.gam;
tstep = 0.05;
t_end = 100;
 
%specify the output points
tspan = [0:tstep:t_end];
p.N=1;
 
[T,S] = ode15s(@(t,x)coscrhs(t,x,p),tspan,init);
C = S(:,1);
Ce = S(:,2);
formatSpecF = '%6.2f\n';
figure(2*j-1) % a phase portrait
plot(c,ce1,'--',c,ce2,'--',C,Ce,C(1),Ce(1),'*')
legend('dc/dt=0','dc_e/dt=0')
xlabel('c')
ylabel('c_e')
axis([0 1.5 0 axlist(j)])
title(strcat('p =',sprintf(formatSpecF,p.p),'\mu M'),'fontsize',18)
ycords=[0.66,0.61;0.8,0.775];
xcords=[0.42,0.52;0.42,0.52];
annotation('textarrow',xcords(j,:),ycords(j,:),'linewidth',2) 
 

figure(2*j) % a time sequence
 plot(T,C,T,Ce,'linewidth',2)
 
xlabel('Time','fontsize',20)
legend('Ca^{++}','Ca^{++}_e')
  title(strcat('p =',sprintf(formatSpecF,p.p),'\mu M'),'fontsize',18)
end


% now integrate the pde
  p.N=600;  % number of spatial grid points
  p.L=30;  
  p.h=p.L/p.N;
 p.p=parlist(2)
  % specify diffusion coefficients
 p.dv = 1;p.Dc; % diffusion coefficient for Ca
 p.du = 0.0; % diffusion coefficient for IP3
 p.sc = [1;2*ones(p.N-2,1);1];
  X = p.h*(1:p.N)';
  %set initial data
   
 V =  exp(-3*X.^2 );  % calcium
 U =p.p*ones(p.N,1);  % IP3: this is in the range where we expect a solitary pulse

 % use the steady state from the ode solution as initial data for the
    % pde
  ce0=Ce(end)*ones(p.N,1);
 

init=[U;V;ce0];
   

 
tstep = 0.05; % time between plots
t_end = 8; %total time to run simulation
 
%specify the output points
tspan = [0:tstep:t_end];

% integrate the system of ode's:
[T,S] = ode23( @(t,x)pdeRHS(t,x,p),tspan,init, odeset('maxstep',1));  

% plot the solution at each time step   
 

for j = fix(length(T)/2):fix(length(T)/2)
   figure(6)
  subplot(2,1,1)
  plot(X,S(j,p.N+1:2*p.N),'linewidth',2)
  xlabel('x','fontsize',20)
   ylabel('Ca^{++}(x,t)','fontsize',20)
  title(strcat('p =',sprintf(formatSpecF,p.p),'\mu M'),'fontsize',18)
  subplot(2,1,2)
  plot(X, S(j,2*p.N+1:3*p.N),'linewidth',2)
   xlabel('x','fontsize',20)
    ylabel('C_e(x,t)','fontsize',20)
     
 figure(7)
 plot(c,ce1,'--',c,ce2,'--',S(j,p.N+1:2*p.N), S(j,2*p.N+1:3*p.N),'linewidth',2)
   xlabel('Ca','fontsize',20)
    ylabel('C_e','fontsize',20)
    legend('dc/dt=0','dc_e/dt=0')
     title(strcat('p =',sprintf(formatSpecF,p.p),'\mu M'),'fontsize',18)
end
 
 
 figure(8)
  
  xlabel('time')
 
  mesh(X,T,S(:,p.N+1:2*p.N))

    % % Now find the speed:
 thresh = 0.5;
 % for each X value find the first time the solution crosses the threshold
 for j = 1:p.N
     jmin = min(find(S(:,p.N+j)>=thresh));
     Tc(j) = T(jmin);
 end
 
  
 q=polyfit(Tc,X,1)
 spest= q(2)+q(1)*Tc;
  figure(9)
 plot(X,Tc,spest,Tc,'--')
 xlabel('x')
 ylabel('t')
  speedest = q(1)

 title(strcat('Speed =',sprintf(formatSpecF,speedest)),'fontsize',18)

%the right hand side for pde (MoL) simulation:
function s_prime=pdeRHS(t,s,p)
 
%There are four variables:  P, C, y, ce 
% only P and C are diffusing
p.h = p.L/p.N;
 
scv = p.dv/p.h^2;
scu = p.du/p.h^2;
P = s(1:p.N);
C = s(p.N+1:2*p.N);
ce = s(2*p.N+1:3*p.N);

% this is method of lines
% evaluate the ode part
%[Fp;Fc;Fy;Fce]=
 
out=coscrhs(t,s(p.N+1:3*p.N),p);
 
Fc=out( 1:p.N);
Fce=out(p.N+1:2*p.N);  % reaction only, no diffusion
FP = scu*(-p.sc.*P+[0;P(1:end-1)]+[P(2:end);0]);  % IP3 is simply diffusing, no reaction term
FC = scv*(-p.sc.*C+[0;C(1:end-1)]+[C(2:end);0]) +Fc ;
 
s_prime = [FP;FC;Fce];
 

%the right hand side for ode simulation:
 
  
function out=coscrhs(t,s,p)
 % evaluate the ode dynamics
% there are two variables
 
c = s(1:p.N); % calcium

ce = s(p.N+1:2*p.N); %ER calcium

Po = p.p*(c.^2./(c.^2+p.ph1^2)).*(p.ph2./(c+p.ph2)); %open probability

Jipr = (p.kf*Po+p.a).*(ce-c);
Jserca = p.ks*c;
Jpm = p.km*c;
Jin = p.a0+p.a1*p.p;

Fc = Jipr-Jserca+p.del*(Jin-Jpm);
 
Fce = p.gam*(Jserca-Jipr);
out = [Fc;Fce];
 
 