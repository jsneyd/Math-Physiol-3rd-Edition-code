% Program to solve the initial stages of the modified Forti et al model of
% phototransduction. For the exercises of Chapter 19, Keener and Sneyd.


clear all
close all
clc

% specify all the parameter values
par.T0 = 1000; par.P0 = 100;
par.alpha = 20; par.beta = 10.6; par.eps = 0.5;  
par.tau1 = 0.1; par.tau2 = 10;

IC = [0 0 0];                   % initial conditions
stim = [1 2 4 100 200 400];     % The stimulus strengths, going up in factors of 2, just for convenience.
tspan=linspace(0,1,500);        % The time interval on which to solve the odes


for i=1:6
[t,U] = ode15s(@(t,y)rhs(t,y,par,stim(i)),tspan,IC);
figure(1);
plot(t,U(:,3),'LineWidth',2')           % plot the responses
hold on
figure(2);
plot(t,U(:,3)/stim(i),'LineWidth',2')   % plot the responses scaled by the stimulus strength.
hold on
end

% make the figures look nicer
figure(1)
set(gca,'FontSize',14)
xlabel('time (s)')
ylabel('P^*')
box off
lgd = legend('1','2','4','100','200','400');
title(lgd,'stimulus strength','FontSize',14)
hold off

figure(2)
set(gca,'FontSize',14)
xlabel('time (s)')
ylabel('P^* (scaled by stimulus strength)')
box off
lgd = legend('1','2','4','100','200','400');
title(lgd,'stimulus strength','FontSize',14)
hold off

% Now save the figures. This bit is just for my convenience. You might want
% to delete these lines

saveas(1,'../../../Math-Physiol-3rd-Edition/figures/chap_19_retina/exercises/visex1_fig3.png')
saveas(2,'../../../Math-Physiol-3rd-Edition/figures/chap_19_retina/exercises/visex1_fig4.png')


%% The ODEs -------------------------------------------------
function dUdt=rhs(t,U,par,stim)

Rs = U(1); Ts=U(2); Ps=U(3);
light = stim*(heaviside(t) - heaviside(t-0.01));  % Stimulus is a square impulse, of width 0.01 and height stim.

dUdt(1) = light - par.alpha*Rs;
dUdt(2) = par.eps*Rs*(par.T0-Ts-Ps) - par.beta*Ts*(Ts+Ps) - par.tau1*Ts*(par.P0-Ps) + par.tau2*Ps;
dUdt(3) = par.tau1*Ts*(par.P0-Ps) - par.tau2*Ps;

dUdt = dUdt';

end




