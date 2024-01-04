%  -------------------------------------------------------------------
%
%   Program to compute the response of the AMPA model to a delta function
%   input
%
%   For Chapter 8, Fig. 8.13 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 

clear all
close all
clc

set(0,                           ...
    'defaultaxesfontsize', 20,   ...
    'defaultaxeslinewidth', 2.0, ...
    'defaultlinelinewidth', 2.0);

%parameters
p.Rb=13.0e6;
p.Ru1=5.9;
p.Ru2=8.6e4;
p.Rd=900.0;
p.Rr=64.0;
p.Ro=2.7e3;
p.Rc=200.0;
p.alpha=1.35;
p.beta=200.0;

tstep = 0.001;
t_end = 0.03;
 
%specify the output pooints
tspan = [0:tstep:t_end];

s0 = [1;0;0;0;0;1];
 
[T,S] = ode15s(@(t,x)rhs(t,x,p),tspan, s0);  
c0 = S(:,1);
c1 = S(:,2);
c2 = S(:,3);
D1 = S(:,4);
D2 = S(:,5);
simpleC = S(:,6);
O = 1-c0-c1-c2-D1-D2;
simpleO = 1-simpleC;
 
plot(1000*T,O,1000*T,simpleO,'--')
xlabel('time(ms)')
ylabel('faction of open channels')
legend('boxoff')
legend('full model','simple model')

%% evaluate the ode dynamics
function out=rhs(t,s,p)

c0 = s(1); 
c1=s(2);
c2=s(3);
D1=s(4);
D2=s(5);
simpleC=s(6);
simpleO=1-simpleC;
O=1-c0-c1-c2-D1-D2;
TT = 1000*((t>0) - (t>0.001));
 
c0p = -p.Rb*TT*c0 + p.Ru1*c1;
c1p = p.Rb*TT*c0 - p.Ru1*c1 + p.Rr*D1 - p.Rd*c1 + p.Ru2*c2 - p.Rb*TT*c1;
c2p = -p.Ru2*c2 + p.Rb*TT*c1 + p.Rr*D2 - p.Rd*c2 + p.Rc*O - p.Ro*c2;
D1p = p.Rd*c1-p.Rr*D1;
D2p = p.Rd*c2 - p.Rr*D2;
simpleCp= -p.alpha*simpleC*TT + p.beta*simpleO;

out = [c0p c1p c2p D1p D2p simpleCp]';
end
