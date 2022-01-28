
% Code to solve a simple version of the Longtin-Milton model of the pupil
% light reflex. This version of the model doesn't have a fixed delay but
% uses instead a four-stage linear filter to create a delay. So it uses
% ode45 to solve, and hence runs fine under Octave also. The parameters of
% the four-stage linear filter are chosen pretty arbitrarily here.

function longtin_milton_ex_nodelay
clear all
close all
clc

IC = [1 1 1 1 1];
tspan = linspace(0,3000,1000);
[t,sol] = ode45(@(t,y)rhs(t,y), tspan,IC);
area = areafun(sol(:,5));
plot(t,sol(:,4),'b-o',t,area,'r-o')
xlabel('Time (s)');
ylabel('Pupil area (mm^2)');
end

%%
function yp = rhs(t,y)
gamma = 8;
I = 10;
phibar = 1;
k1 = 0.1; k2=1;

x = y(5);
x1 = y(1);
x2 = y(2);
x3 = y(3);
x4 = y(4);

yp(1) = k1*(x - k2*x1);
yp(2) = k1*(x1 - k2*x2);
yp(3) = k1*(x2 - k2*x3);
yp(4) = k1*(x3 - k2*x4);
yp(5) = gamma*capF(log(I*areafun(x4)/phibar))-x;

yp = yp';
end


%% 
function out=capF(x)
out = heaviside(x)*x;
end

%%
function out = areafun(x)
lam = 30;
n = 7;   % vary n to get different oscillatory behaviour
theta = 10;
out = lam*theta^n./(x.^n + theta^n);
end
