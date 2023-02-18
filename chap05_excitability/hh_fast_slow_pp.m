
function hh_fast_slow
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 


global Iapp gna gk gl vna vk vl  
global dumv

Iapp=40.0;
gna=120.0; gk=36.0; gl=0.3;
vna=115.0; vk=-12.0; vl=10.6; cm=1.0;

vlist = linspace(vk+0.01,vna-0.01,100)
[n0 nbar]=get_nullclines(vlist);
figure(1)

plot(vlist,n0,vlist,nbar)
xlabel('v')
ylabel('n')
axis([-20 120 0 1])
hold on

output=[vlist' n0' nbar'];
save fastslow.dat output -ascii
end

%%
function [n ,nbar]=get_nullclines(vlist)
global Iapp gna gk gl vna vk vl  
global dumv
n0=0;  % initial guess for the solver

for i=1:length(vlist)
    dumv=vlist(i);
    
    n(i)=fminsearch(@nullcline,n0); % this uses a Matlab routine to find the minimum of a function
    n0=n(i);
    alphan=0.01*(10.0-dumv)./(exp((10.0-dumv)/10.0)-1.0);
    betan=0.125*exp(-(dumv)/80);
    nbar(i)=alphan/(alphan+betan);
end

end

%%
function out=nullcline(n)
global Iapp gna gk gl vna vk vl  
global dumv

alpham=0.1*(25.0-dumv)/(exp((25.0-dumv)/10.0)-1.0);
betam=4.0*exp(-(dumv)/18.0);
m=alpham/(alpham + betam);

out=((n^4)*gk*(dumv-vk) + gna*(m^3)*(0.8-n)*(dumv-vna) + gl*(dumv-vl) - Iapp)^2;
end