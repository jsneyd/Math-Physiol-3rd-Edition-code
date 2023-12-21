
clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0,...
   'defaultlinemarkersize',12);

global p
global O U Otilde Utilde   % keep them to use as the starting point for the nonlinear solve
global ic

%%%
% starting initial conditions for the nonlinear solve. This is needed   
% as otherwise the nonlinear solve can get weird  (confused), if you start at a
% bad place. It's highly inefficient to keep these as global variables, but
% oh well, it seems to work.

O=0; U=0; Otilde=0; Utilde=0;   
%%%
% normalize time scale
tscale = 249/120;

 formatSpecF = '%6.2f\n';
p.sigma = 1.4;
p.sigmaCO = 33;
p.beta = 200;
p.KO = 30*p.sigma;
p.THb = 2200;
p.RT = 1.7e-2;   % units of uM
p.D = 32.5;  % same for O2 and CO
 
p.PO = 0.21*(760-47);  % the 47 is from the water vapor
 
% with these parameters, O2 = 104/40 mm Hg with no CO

p.D = 32.5;
 
p.gamma = 37;

PO_newList = [150,1*(760) - 47, 2.5*(760) - 47];
 
for  ic = 1:3
p.ic = ic;
p.PO_new= PO_newList(ic);
p.PO = 150;  % initial O2 partial pressure
p.PCO = 1;  % initial partial pressure of CO
% initialize the nonlinear solve
O=p.sigma*40;
 U=0; 
 Otilde=104*p.sigma ;
 Utilde=0; 
 
% solve the ODE for TO
%init = [5000 5000];
init = [6800,0];
tspan = linspace(0,420,1000);
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[T,Y] = ode15s(@Tode,tspan,init,options);
TO=Y(:,1);
TU=Y(:,2);
for j = 1:length(T)
    getO_Otilde_U_Utilde(TO(j),TU(j));
    svO(j) = O;
    svU(j) = U;
    svOtilde(j) = Otilde;
    svUtilde(j) = Utilde;
end

figure( ic)
plot(T*tscale,TO,'r',T*tscale,TU,'g')
xlabel('time (minutes)')
ylabel('total concentration (\mu M)')
 
box off
legend('O_2','CO')
if (ic==1)
  
  title( 'P_{O_2} = 0.21 Atm')
end
if (ic ==2)
     title( 'P_{O_2} = 1 Atm')
end
if (ic ==3)
     title( 'P_{O_2} = 2.5 Atm')
end
% Now find the point where CO is half of its maximum (as it decays)
[M,J] = max(TU);
half = find(TU > M/2,1,'last');
half_clearance = (tspan(half)-20) *tscale
end

%% ODE for TO
function out = Tode(t,x)
TO = x(1);
TU = x(2);
global p 
global O U Otilde Utilde   % these are used as the starting point for the nonlinear solve

if (t>20)   % put into a zero CO and high O2 environment after t = 20
    p.PCO = 0;
    p.PO = p.PO_new;

end
getO_Otilde_U_Utilde(TO,TU);
out(1) = p.D*(p.sigma*p.PO - Otilde) - p.gamma*O;
out(2) = p.D*(p.sigmaCO*p.PCO - Utilde);
out = out';

end

%% get O, Otilde, U, Utilde as functions of TO and TU

function getO_Otilde_U_Utilde(TO,TU)
global p
global O U Otilde Utilde   % keep them to use as initial conditions

options = optimset('Display','off');
% first find O and U from TO and TU
solOU = fsolve(@(x)getOU(x,TO,TU),[O U],options);  % starting point is the solution from the last solve.
O = solOU(1);
U = solOU(2);

soltildes = fsolve(@(x)getOtildeUtilde(x,TO,TU),[Otilde Utilde],options);
Otilde = soltildes(1);
Utilde = soltildes(2);

end

%%
function out = getOU(x,TO,TU)
O = x(1);
U = x(2);
out(1) = concO(O,U) - TO;
out(2) = concU(O,U) - TU;

end

%%
function out = getOtildeUtilde(x,TO,TU)
global p

Otilde = x(1);
Utilde = x(2);
out(1) = p.D *(p.sigma*p.PO - Otilde) + TO - concO(Otilde,Utilde);
out(2) = p.D *(p.sigmaCO*p.PCO - Utilde) + TU - concU(Otilde,Utilde);

end

%% oxygen saturation function
function out=concO(O,U)
% total O as function of O and U
global p

out = O + 4*p.THb*O.^4./(p.KO^4 + O.^4 + p.beta^4*U.^4);
end

%% CO saturation function
% total U as function of O and U
function out=concU(O,U)
global p

out = U + 4*p.THb*p.beta^4*U.^4./(p.KO^4 + O.^4 + p.beta^4*U.^4);
end
