
function main

global Iapp gna gk gl vna vk vl cm
global dumv

Iapp=50.0;
gna=120.0; gk=36.0; gl=0.3;
vna=115.0; vk=-12.0; vl=10.6; cm=1.0;

n0=0;
for i=1:200
    dumv=-12+127*i/200;
    v(i)=dumv;
    n(i)=fminsearch(@nullcline,n0);
    n0=n(i);
    alphan=0.01*(10.0-dumv)./(exp((10.0-dumv)/10.0)-1.0);
    betan=0.125*exp(-(dumv)/80.0);
    nbar(i)=alphan/(alphan+betan);
end
plot(v,n,v,nbar)
output=[v' n' nbar'];

save fastslow.dat output -ascii

%%%%%%%%%%%%%%%%%%%%%%%%
function out=nullcline(n)
global Iapp gna gk gl vna vk vl cm
global dumv

alpham=0.1*(25.0-dumv)/(exp((25.0-dumv)/10.0)-1.0);
betam=4.0*exp(-(dumv)/18.0);
m=alpham/(alpham + betam);

out=((n^4)*gk*(dumv-vk) + gna*(m^3)*(0.8-n)*(dumv-vna) + gl*(dumv-vl) - Iapp)^2;