close all
clear all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global tau


% First plot a typical time series
p = 0.5;
tau = 1;
Ninit = 1;
N = Ninit;
T = 0;
delt = 1;

for count = 1:7
    tt = linspace((count-1)*delt,count*delt,50);
    [Ti,Ni] = ode15s(@(t,n)rhs(t,n),tt,Ninit);
    Ninit = Ni(end)*(1-p);
    N = [N,Ni'];
    T = [T,Ti'];
end
figure(1)
plot(T,N,'LineWidth',2)
xlabel('t(s)')
ylabel('n')

% now plot the frequency versus mean response
fcount = -2:0.1:2;
for i = 1:size(fcount')
    freq(i) = 10^fcount(i);
    delt = 1/freq(i);
    nbar(i) = (delt - p*tau + (p^2*tau)/(p + exp(delt/tau) - 1))/delt;  % average value
    testfit(i) = 1/(1+freq(i)*tau*p);
end
figure(2)
semilogx(freq,nbar,freq,testfit,'LineWidth',2)
ylim([0,1])


save('STSD.mat','T','N','freq','nbar')


%%
function out = rhs(t,n)
global tau

out = (1-n)/tau;

end
