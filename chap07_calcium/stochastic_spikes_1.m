
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

global x0
N = 50000;
delt = 0.01;
k = 0.5;
x0 = 0.01;  % The system can be reactivated only when x<x0

nthresh = 5; %# of different runs % use 40 for a mu/sigma plot, 
% but it takes a lot of time
thresh = linspace(0.997,0.9995,nthresh); 
% for a mu/sigma plot
% thresh = linspace(0.999,0.999,nthresh);           % for a single run
 mu = zeros(nthresh,1);                              % initialization
sig = mu;                                           % initialization

for j = 1:nthresh
    threshold = thresh(j);
    In = rand(N,1);
    init = 0;  % first starting point
    x1=0;
   
    Tkeep = [0];
    i = 0;
T = [0];
tend = T(end);
Xkeep = [0];

    while i<N
        % find the next spike time
        nn = find(In(i+1:end)>threshold);
        if(~isempty(nn))
        ndx = min(nn);  % this is how many time steps to take
newt = [1:ndx]*delt;

% integrate the de until the next spike time:
 % tspan = [0 newt];
        %[dT,X] = ode15s(@(t,x)rhs(t,x,k),tspan,x1);
Xkeep  = [Xkeep x1*exp(-k*newt)]; %no need to integrate when the analytical solution is known
T = [T tend+newt];
tend = T(end);
Tkeep = [Tkeep tend];
x1=Xkeep(end)+1;  % add the jump
 %here we are using that decay is exponential.  For a more complicated de,
 %numerical integration is required
 % integrate the de until threshold x0 is reached
 %tspan = [0:delt:100]; % we only use a fraction of these times
  %[dT,X] = ode15s(@(t,x)rhs(t,x,k),tspan,x1, options);

 dt=-log(x0/x1)/k; %time to get back to recovery, since the analytical solution is known
 
 ndt = fix(dt/delt); % number of time step to get back to recovery
 newt = [1:ndt]*delt;
 Xkeep = [Xkeep x1*exp(-k*newt)];
 T = [T tend+newt];
 tend = T(end); % x will be below threshold on the next time step
 i = length(T);
 x1 = Xkeep(end);  % this is the starting value for the next integration
        else
            i=N;
        end
    end
 
   
    figure(j)
    plot(T,Xkeep)
    xlabel('time')
    ylabel('x')
    ISI = Tkeep(2:end) - Tkeep(1:end-1);            % the interspike intervals
    mu(j) = mean(ISI);
    sig(j) = std(ISI);
end

figure(10)
plot(mu,sig,'o')
hold on

% Add the linear fit to the plot
p = polyfit(mu,sig,1)
xp = linspace(9,35,100);
plot(xp,polyval(p,xp));
xlabel('\mu')
ylabel('\sigma')

%% ode equations

function out = rhs(t,x,k)
out = -k*x;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f,direction, isterminal] = event(t,x)

global x0
direction = 1;
f = x-x0;
isterminal = 1;
end
