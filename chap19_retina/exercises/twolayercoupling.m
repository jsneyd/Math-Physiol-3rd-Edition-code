
% Program to find the response to a Heaviside function, when there is
% current spread in both photoreceptor and horizontal cell layers.

clear all 
close all
clc

lh = 0.05;   % voltage spread in horizontal cell layer 
lr = 0.05;   % voltage spread in receptor cell layer
k = 10;
L = 1;      % Solve from -L to L

% Find delta and gamma

a = 1+k;
b = -lh^2-lr^2;
c = lh^2*lr^2;

delsq = (-b + (b^2-4*a*c)^0.5)/(2*a);
delta = delsq^0.5;

gamsq = lr^2*lh^2/((1+k)*delsq);
gamma = gamsq^0.5;

% Now solve in two bits, positive x and negative x.
x1 = linspace(0,L,100);
Wgam = gamma^2*(1 - 0.5*exp(-x1/gamma));
Wdel = delta^2*(1 - 0.5*exp(-x1/delta));
V1 = (Wdel - Wgam)/((1+k)*(delta^2-gamma^2));

x2 = linspace(-L,0,100);
Wgam = gamma^2*exp(x2/gamma)/2;
Wdel = delta^2*exp(x2/delta)/2;
V2 = (Wdel - Wgam)/((1+k)*(delta^2-gamma^2));

plot(x1,V1,'r',x2,V2,'r','LineWidth',2)
xlabel('x')
ylabel('V')
set(gca,'FontSize',14)
