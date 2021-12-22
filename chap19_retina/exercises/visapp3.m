
clear all
close all
clc

% set parameters
f0 = 1;
epsilon = 0.2;
w = 1;
a = 2;

% numerical solve of ode. Note we use i*w*t but Matlab only plots the real
% parts, which is all we are interested in
tspan = [0 25];
xinit = 1;
[t,x] = ode45(@(t,x) 1 + epsilon*exp(i*w*t)-a*x*x, tspan, xinit);

% Now calculate the first-order frequency response
x0 = sqrt(f0/a);
x1 = epsilon/(2*sqrt(a*f0) + i*w);
fr = x0 + x1*exp(i*w*t);

% Plot them both in the same graph
plot(t,x,t,fr,'Linewidth',2)
legend('exact','first order frequency response')
box off
ylabel('x')
xlabel('t')
set(gca,'Fontsize',14)
title(['\epsilon = ',num2str(epsilon)])