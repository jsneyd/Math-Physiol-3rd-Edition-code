% Reichardt motion detector, responding to a square wave pulse.

clear all
close all
clc

dx = 3;
dt = 6;
t = linspace(-15,15,200);

%% First do the single component
figure(1)
plot(t,s(0,t))
hold on
plot(t,s(dx,t))
plot(t,s(0,t-dt))
plot(t,s(dx,t).*s(0,t-dt),'--','LineWidth',2)
hold off

%% Next do the full detector

figure(2)
plot(t,s(0,t),t,s(dx,t-dt))
hold on
plot(t,s(0,t-dt),t,s(dx,t))
figure(3)
plot(t,s(0,t).*s(dx,t-dt),t,s(0,t-dt).*s(dx,t))
plot(t,s(0,t).*s(dx,t-dt)-s(0,t-dt).*s(dx,t),'LineWidth',2)


%%
function out = s(x,t)
c=0.5;
s0=0;
    out = s0 + exp(-(x-c*t).^2);
end