function ex11_fig1
clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

%%
% First do the usual QSSA solution, with e0/s0 small

% Set up the parameters
par.e0 = 0.1;
par.s0 = 1;
par.k1 = 1;
par.km1 = 1;
par.k2 = 1;
par.Km = (par.km1+par.k2)/par.k1;

eps = par.e0/(par.s0)                   % Display to screen, just for fun

tspan = linspace(0,30,5000);            % time interval for the solution
IC = [par.s0 0];                        % initial condition
[t,U] = ode15s(@(t,y)rhs(t,y,par),tspan,IC);
figure(1)
plot(U(:,1),U(:,2),'LineWidth',2)       % plot the phase plane
hold on

% Now calculate the slow manifold and add it to the plot
s = linspace(0,1.2,200);
c = par.e0*s./(par.Km+s);
plot(s,c,'r--','LineWidth',2)

% Make the plot look pretty
set(gca,'FontSize',16)
xlabel('substrate concentration (any units you like)')
ylabel('complex concentration (any units you like)')
legend('solution','slow manifold','Location','NorthWest')
box off

%%
% Next do the modified QSSA solution, with e0/(s0+Km) small

% Set up the parameters again, using different values to make the new epsilon
% small, but with e0/s0 not small.
par.e0 = 1;
par.s0 = 1;
par.k1 = 1;
par.km1 = 20;
par.k2 = 20;
par.Km = (par.km1+par.k2)/par.k1;

eps = par.e0/(par.s0+par.Km)            % Display to screen, just for fun.

tspan = linspace(0,30,5000);            % time interval for the solution
IC = [par.s0 0];                        % initial condition
[t,U] = ode15s(@(t,y)rhs(t,y,par),tspan,IC);
figure(2)
plot(U(:,1),U(:,2),'LineWidth',2)       % plot the phase plane
hold on

% Now calculate the slow manifold and add it to the plot
s = linspace(0,1.2,200);
c = par.e0*s./(par.Km+s);
plot(s,c,'r--','LineWidth',2)

% Make the plot look pretty
set(gca,'FontSize',16)
xlabel('substrate concentration (any units you like)')
ylabel('complex concentration (any units you like)')
legend('solution','slow manifold','Location','NorthWest')
box off

% Now save the figures. This bit is for my convenience. You might want
% to delete these lines

saveas(1,'../../../Math-Physiol-3rd-Edition/figures/chap_01_biochem/exercises/ex11_fig1.png')
saveas(2,'../../../Math-Physiol-3rd-Edition/figures/chap_01_biochem/exercises/ex11_fig2.png')


%% ODEs

function out=rhs(t,y,par)

s = y(1);
c = y(2);

out(1) = par.km1*c - par.k1*s*(par.e0-c);
out(2) = par.k1*s*(par.e0-c) - (par.k2+par.km1)*c;

out = out';

 