% code to simulate respiration control (Chapter 14)


function respiration_control
 
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
    
 %parameters for repiratory control
p.tau = 1.0;
p.m = 0.5;
p.mu = 5.0;
p.k = 1.0;
p.E1 = 0.5;
p.E2 = 0.3;
p.K = 0.2;
 
 

% integrate the ode
 init = [0.2,0.2,0,0.]; %initial data for the odes
tstep = 0.1;
t_end = 100;
 
%specify the output points
tspan = [0:tstep:t_end];
 
[T,S] = ode15s(@(t,x)rhs(t,x,p),tspan,init);
I1 = S(:,1);
I2 = S(:,2);
x = S(:,3);
xdot = S(:,4);
 
 
figure(1)
 plot(T,I1,T,I2,'--','linewidth',2)
 legend('I_1','I_2')
xlabel('Time','fontsize',20)
 
figure(3)
plot(T,x)
%the right hand side for ode simulation:
   
function out=rhs(t,s,p)
 % evaluate the ode dynamics
  
I1 = s(1); %
I2 = s(2); %  
x = s(3);  %  
xdot = s(4);  %

arg1=p.E1-I2;
F1 = (arg1>0)*2*arg1^2/(p.K+arg1);

f= 2.5*x^3/(1+x^3);
 
arg2 = p.E2-I1+f;
F2 = (arg2>0)*2*arg2^2/(p.K+arg2);
 
FI1 =  (F1-I1)/p.tau;
FI2 = (F2-I2)/p.tau;
Fx = xdot;
Fxd = (I1-p.k*x-p.mu*xdot)/p.m;
 
out = [FI1;FI2;Fx;Fxd];
 