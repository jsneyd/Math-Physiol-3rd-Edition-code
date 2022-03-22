function main

close all
clear all
clc

n = 200 ; % number of coupled oscillators

% First plot the uncoupled case, to see the natural frequencies
G=0;
tspan = linspace(0,60,1200);
x0(1:n) = 0.1; %linspace(0,0.8,n);
x0(n+1:2*n) = 0;
[t,sol] = ode15s(@(t,x) FHN(t,x,n,G),tspan,x0);
%plot(t,sol(:,2),t,sol(:,100),t,sol(:,200))

% Find the frequency of each oscillator
for i=1:n
    [~,idx] = findpeaks(sol(:,i));
    period(i) = (t(idx(end)) - t(idx(end-10)))/10; % average over 10 cycles to get smoother result
    frequency(i) = 1/period(i);
end
plot(frequency)
hold on


% Now plot the coupled case
G=50;
[t,sol] = ode15s(@(t,x) FHN(t,x,n,G),tspan,x0);

% Find the frequency of each oscillator
for i=1:n
    [~,idx] = findpeaks(sol(:,i));
    period(i) = (t(idx(end)) - t(idx(end-10)))/10; % average over 10 cycles to get smoother result
    frequency(i) = 1/period(i);
end
plot(frequency)

end

% -------------------------------------------------------------------------------
%% FHN equations
function dxdt = FHN(t,x,n,G)

v = x(1:n);
w = x(n+1:2*n);

alpha = 0.4;
eps = linspace(0.02,0.08,n);

dxdt(n+1:2*n) = v - alpha';  % no coupling in w
for i=2:n-1
    dxdt(i) = (1/eps(i))*(v(i).*(v(i)-1).*(alpha - v(i)) - w(i)) + G*(v(i+1)-v(i)) + G*(v(i-1)-v(i));
end
dxdt(1) = (1/eps(i))*(v(1).*(v(1)-1).*(alpha - v(1)) - w(1)) + G*(v(2)-v(1));
dxdt(n) = (1/eps(i))*(v(n).*(v(n)-1).*(alpha - v(n)) - w(n)) + G*(v(n-1)-v(n));

dxdt = dxdt';
end