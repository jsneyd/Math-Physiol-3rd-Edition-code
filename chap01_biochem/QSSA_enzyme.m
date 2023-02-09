function enzyme
% asymptotic analysis of enzyme kinetics equations - QSSA
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
   
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
xlabel('\eta','fontsize',20)
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

 

function out = rhsqss(t,y,par)
s = y(1);
x = y(2);
out(1) = -s+x*(s+par.a) ;
out(2) = (s-x*(s+par.k))/par.e;
out = out';  % don't forget this or you get rubbish.
end

