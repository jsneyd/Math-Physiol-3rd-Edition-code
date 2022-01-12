% Reichardt motion detector, responding to a square wave pulse.

clear all
close all
clc

dx = 2.8;
dt = 6;
t = linspace(-15,15,200);
c = -0.5;    % preferred direction
s0 = 1;     % not a raised stimulus

%% First do the single component

figure(1)
plot(t,s(0,t,c,s0),'LineWidth',2)
hold on
plot(t,s(dx,t,c,s0),t,s(0,t-dt,c,s0),'LineWidth',2)
plot(t,s(dx,t,c,s0).*s(0,t-dt,c,s0),'--','LineWidth',2)
hold off
xlabel('t')
ylabel('response')
box off
set(gca,'FontSize',14)
if (c>0) 
    title('preferred direction') 
else
    title('null direction')
end
legend('stimulus','right-hand receptor','left-hand receptor','combined response','Location','northwest')


%% Next do the full detector

figure(2)
plot(t,s(0,t,c,s0).*s(dx,t-dt,c,s0),t,s(0,t-dt,c,s0).*s(dx,t,c,s0),'LineWidth',2)
hold on
plot(t,s(0,t,c,s0).*s(dx,t-dt,c,s0)-s(0,t-dt,c,s0).*s(dx,t,c,s0),'LineWidth',2)
xlabel('t')
ylabel('response')
box off
set(gca,'FontSize',14)
if (c>0) 
    title('preferred direction') 
else
    title('null direction')
end
legend('first product', 'second product','combined response','Location','northwest')


%%
function out = s(x,t,c,s0)
    out = s0 + exp(-(x-c*t).^2);
end