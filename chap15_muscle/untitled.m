
clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0,...
   'defaultlinemarkersize',12);

w0 = linspace(0,10,100);
D1 = 1;
Dc_over_D1 = 0.1/(1+0.1);
v1 = 2*Dc_over_D1*(exp(w0)-1)./(exp(w0)+1);

plot(w0,v1,'--r')
hold on

f = @(w0,v) getf(w0,v);
fcontour(f, [0 10 0 0.4],'-b','LevelList',[0],'LineWidth',2);

xlabel('\omega_0')
ylabel('v\delta/D_1')

%%
function out = getf(w0,v)
    wl = 10*v;
    out = v - wl.^2.*(1 - exp(w0-wl))./( (exp(w0)-1).*(exp(-wl)-1) + wl.*(exp(w0-wl)-1) );
end