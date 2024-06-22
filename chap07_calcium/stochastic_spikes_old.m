
%  -------------------------------------------------------------------
%
%   Compute a stochastic train of spikes in a simple model.
%
%   For Chapter 7, Section 7.10.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 

clear all
close all
clc
set(0,                           ...
    'defaultaxesfontsize', 20,   ...
    'defaultaxeslinewidth', 2.0, ...
    'defaultlinelinewidth', 2.0);

N = 50000;
delt = 0.01;
k = 0.5;
x0 = 0.01;  % The system can be reactivated only when x<x0

nthresh = 1; %# of different runs % use 40 for a mu/sigma plot, 
% but it takes a lot of time
thresh = linspace(0.997,0.9995,nthresh);          % for a mu/sigma plot
% thresh = linspace(0.999,0.999,nthresh);           % for a single run

mu = zeros(nthresh,1);                              % initialization
sig = mu;                                           % initialization

for j = 1:nthresh
    threshold = thresh(j);
    In = rand(N,1);
    init = 0;  % first starting point
    Tkeep = 0;
    
    for i = 1:N
        tspan = linspace(0,delt,2);
        [T,X] = ode15s(@(t,x)rhs(t,x,k),tspan,init);
        init = X(end);
        if (In(i)>threshold & init<x0)
            init = X(end)+1;
            Tkeep = [Tkeep i*delt];                 % record the time of the spike
        end
        Xkeep(i) = init;
    end
    figure(j)
    plot(Xkeep)
    xlabel('time')
    ylabel('x')
    ISI = Tkeep(2:end) - Tkeep(1:end-1);            % the interspike intervals
    mu(j) = mean(ISI)
    sig(j) = std(ISI)
end

figure(10)
plot(mu,sig,'o')
hold on

% Add the linear fit to the plot
p = polyfit(mu,sig,1)
xp = linspace(9,35,100);
plot(xp,polyval(p,xp));


%% ode equations

function out = rhs(t,x,k)
out = -k*x;
end