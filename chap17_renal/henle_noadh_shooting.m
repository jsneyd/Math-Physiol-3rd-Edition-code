% this code is a so far a failed attempt to solve the no ADH case using shooting.
clear all
close all
clc
 
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global p

p.P = 0.9;
p.DPd = 0.15;
p.Hd = 0.1;
p.rhod = 0.15;

init = 1;
Qa = -0.2027;
 
Qa = -1.3
%yspan = linspace(0,1,1000);
%[Y,Sdsol] = ode15s(@(y,Sd)rhs(y,Sd,Qa),yspan,init);
yspan = [0:.1:0.1];
for j = 1:length(yspan)
    Y = yspan(j)

Sdsol = [0:.002:2];
out(:,j) = rhs(Y,Sdsol,Qa);

 figure(j)
 plot(Sdsol,out(:,j))

 end
 
size(out)

    tm = (1-Sdsol-p.DPd*p.Hd*Y)/(p.rhod*p.Hd);
    Qd = 1 + tm;
    Qs = -Qd-Qa;
    Cd = Sdsol/Qd;
    Ca=Cd-p.P*(Y-1)/Qa;
    Sa =Qa.*Ca;
    Ss = -Sa - Sdsol;
    Qd = 1 + (1-Sdsol - p.DPd*p.Hd*Y)/(p.rhod*p.Hd);
    
%plot(  Y, Qs,Y, Qd)
%test = Sdsol(end) - (1 - p.DPd*p.Hd + p.rhod*p.Hd*(Qa+1))
box off
%% the ODE
function out = rhs(y,Sd,Qa)
    global p

    %Sd1 = 1 - p.DPd*p.Hd + p.rhod*p.Hd*(Qa+1);

    tm = (1-Sd-p.DPd*p.Hd*y)/(p.rhod*p.Hd);
    Qd = 1 + tm;
    Qs = -Qd-Qa;
    Cd = Sd./Qd;
    Ca=Cd-p.P*(y-1)/Qa;
    Sa =Qa*Ca;
    Ss = -Sa - Sd;
    Qd = 1 + (1-Sd - p.DPd*p.Hd*y)/(p.rhod*p.Hd);

    out = p.Hd*(Ss./Qs - Sd./Qd);
end

