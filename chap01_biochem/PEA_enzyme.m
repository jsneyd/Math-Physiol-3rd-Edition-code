function enzyme
% asymptotic analysis of enzyme kinetics equations
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

% Part 1:  Equilibrium approximation
par.a = 2.2;
par.b = 1.7;
par.e = 0.1;

tspan = linspace(0,20,2000);
initial = [1,1];
[t,Y] = ode15s(@(t,y)rhseq(t,y,par),tspan,initial);
s=Y(:,1);
z=Y(:,2);
x = (z-s)/par.a;   % the original variable, x
figure(1)
plot(t,s,t,z,'linewidth',2)
legend('\sigma','z','fontsize',18)
xlabel('\tau','fontsize',20)
axis([0 2 0 1])
% now plot the slow manifold and the solution together
slow_s = linspace(0,1,100);
slow_z = slow_s + par.a*par.b*slow_s./(1+par.b*slow_s);
slow_x = (slow_z-slow_s)/par.a;
figure(2)
plot(s,z,slow_s,slow_z,'--','linewidth',2)
xlabel('\sigma','fontsize',20)
ylabel('z','fontsize',20)
figure(3)
plot(s,x,slow_s,slow_x,'--','linewidth',2)
xlabel('\sigma','fontsize',20)

ylabel('x','fontsize',20)
forplotting1 = [t s z x];
save('temp1.dat','forplotting1')
forplotting2 = [slow_s ; slow_z;  slow_x]';
save('temp2.dat','forplotting2')


% Part 2:  QSS approximation
par.a = 0.5;
par.k = 1.5;
par.e = 0.05;

tspan = linspace(0,20,2000);
initial = [1,0];
[t,Y] = ode15s(@(t,y)rhsqss(t,y,par),tspan,initial);
s=Y(:,1);
x=Y(:,2);
 
figure(4)
plot(t,s,t,x,'linewidth',2)
legend('\sigma','z','fontsize',18)
xlabel('\tau','fontsize',20)
axis([0 4 0 1])
% now plot the slow manifold and the solution together
slow_s = linspace(0,1.2,100);
slow_x = slow_s./(slow_s+par.k) ;
figure(5)
plot(s,x,slow_s,slow_x,'--','linewidth',2)
xlabel('\sigma','fontsize',20)
ylabel('x','fontsize',20)
axis([0 1.2 0 0.5])
 
forplotting3 = [t s x];
save('temp3.dat','forplotting3')
forplotting4 = [slow_s ;  slow_x]';
save('temp4.dat','forplotting4')


end


function out = rhseq(t,y,par)
s = y(1);
z = y(2);
out(1) = z-s - par.b*s*(par.a + s - z);
out(2) = par.e*(s-z);
out = out';  % don't forget this or you get rubbish.
end

function out = rhsqss(t,y,par)
s = y(1);
x = y(2);
out(1) = -s+x*(s+par.a) ;
out(2) = (s-x*(s+par.k))/par.e;
out = out';  % don't forget this or you get rubbish.
end

