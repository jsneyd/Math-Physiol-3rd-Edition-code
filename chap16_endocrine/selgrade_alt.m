
%    -------------------------------------------------------------------
%
%     Solution of the Selgrade menstrual cycle model.
%
%     For Chapter 16, Section 16.3.1 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function selgrade
% this uses the Matlab routine dde23 to solve delay differential equations
close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

% parameters
 
par.kL = 2.49;
par.cL = 14;
par.V0L = 1263;
par.V1L = 91000;
par.KmL = 360;
par.KiL = 31.22;
par.cLE = 0.0049;
par.cLP = 0.07;
par.tau = 2;

par.VF=5700.0;
par.cF=8.2;
par.kF=7.3;
par.cFE=0.16;
par.KiF=641.0;
par.cFP=644.0;
par.v = 2.5;
par.b=0.004;
par.c1=0.006;
par.c2=0.05;
par.c3=0.004;
par.c4=0.006;
par.c5=1.3;
par.k=1.43;
par.alpha = 0.77;
par.beta = 0.16;

par.e0 = 48;
par.e1 = 0.1;
par.e2=0.17;
par.e3=0.23;
par.p=0.05;
par.h0=274;
par.h1=0.5;
par.h2=0.5;
par.h3=2;
 

lags = par.tau; 
tspan = linspace(0,100,1000);
sol = dde23(@(t,u,z)ddeRHS(t,u,z,par),lags, @history, tspan) ;
LH = sol.y(2,:); 
SeF = sol.y(6,:);
L4 = sol.y(13,:);
PrF = sol.y(7,:);
 
figure(1)
plot(sol.x,LH)
xlabel('t (days)')
ylabel('LH (\mu g/L)')
box off

Est=par.e0 + par.e1*SeF + par.e2*PrF + par.e3*L4;
figure(2)
plot(sol.x,Est)
xlabel('t (days)')
ylabel('Estradol  (ng/L)')
box off 
end
 
%% RHS of ode  
function out = ddeRHS(t,u,z,par)
LHR =u(1);
LH = u(2);
FSHR = u(3);
FSH = u(4);
RcF = u(5);
SeF = u(6);
PrF = u(7);
S1 = u(8);
S2 = u(9);
L1 = u(10);
L2 = u(11);
L3 = u(12);
L4 =u(13);
L3d = z(12);
L4d = z(13);
PrFd = z(7);


E  = par.e0 + par.e1*SeF + par.e2*PrF + par.e3*L4;
PP  = par.p*(L3 + L4);
PPd = par.p*(L3d + L4d);
Id  = par.h0 + par.h1*PrFd + par.h2*L3d + par.h3*L4d;


LHRp = (par.V0L + par.V1L*(E^8)/(par.KmL^8+E^8))/(1.0+PPd/par.KiL) ...
    - par.kL*(1.0+par.cLP*PP)*LHR/(1.0+par.cLE*E );
LHp = (1.0/par.v)*par.kL*(1.0+par.cLP*PP )*LHR/(1.0+par.cLE*E ) -par.cL*LH;
FSHRp = par.VF/(1.0+Id/par.KiF) - par.kF*(1.0+par.cFP*PP )*FSHR/(1.0+par.cFE*(E^2));
FSHp = (1.0/par.v)*par.kF*(1+par.cFP*PP )*FSHR/(1.0+par.cFE*(E^2)) - par.cF*FSH;
RcFp = par.b*FSH + (par.c1*FSH - par.c2*(LH^par.alpha))*RcF;
SeFp = par.c2*(LH^par.alpha)*RcF + (par.c3*(LH^par.beta) - par.c4*LH)*SeF;
PrFp = par.c4*LH*SeF - par.c5*PrF;
S1p = par.c5*PrF - S1/par.k;
S2p = (S1 - S2)/par.k;
L1p = (S2 - L1)/par.k;
L2p = (L1 - L2)/par.k;
L3p = (L2 - L3)/par.k;
L4p = (L3 - L4)/par.k;
out = [LHRp;LHp;FSHRp;FSHp;RcFp;SeFp;PrFp;S1p;S2p;L1p;L2p;L3p;L4p];

end

function s = history(t)
 
 LHR0=600
LH0=25.34
FSHR0=352
FSH0=142.5
RcF0=9
SeF0=1
PrF0=1
S10=1
S20=1
L10=1
L20=1
L30=1
L40=1
u0 = [ LHR0,LH0,FSHR0,FSH0,RcF0,SeF0,PrF0,S10,S20,L10,L20,L30,L40];
  s = u0';
end