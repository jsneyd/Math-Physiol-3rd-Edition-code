function enzyme

clear all
close all
clc

par.a = 2.2;
par.b = 1.7;
par.e = 0.1;

tspan = linspace(0,20,2000);
initial = [1,1];
[t,Y] = ode15s(@(t,y)rhs(t,y,par),tspan,initial);
s=Y(:,1);
z=Y(:,2);
x = (z-s)/par.a;   % the original variable, x
plot(t,s,t,z)

% now plot the slow manifold and the solution together
slow_s = linspace(0,1,100);
slow_z = slow_s + par.a*par.b*slow_s./(1+par.b*slow_s);
figure(2)
plot(s,z,slow_s,slow_z)

% now plot the slow manifold in the s-x plane
slow_x = par.b*slow_s./(1+par.b*slow_s);
figure(3)
plot(s,x,slow_s,slow_x)

forplotting1 = [t s z x];
save('temp1.dat','forplotting1')
forplotting2 = [slow_s ; slow_z ; slow_x]';
save('temp2.dat','forplotting2')

end


function out = rhs(t,y,par)
s = y(1);
z = y(2);
out(1) = z-s - par.b*s*(par.a + s - z);
out(2) = par.e*(s-z);
out = out';  % don't forget this or you get rubbish.
end
