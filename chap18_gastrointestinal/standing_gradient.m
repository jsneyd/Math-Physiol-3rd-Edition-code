% code to solve the standing gradient water transport equations via
% shooting and bisection

function bisect_function
global D P r c0 alp L N0
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

% guess a lower bound for the unknown initial value
a = 0.1;

% guess an upper bound for the unknown initial value
b = 2;

%call the function bisect(a,b,@fval)
root = bisect(a,b,@transport_eqs)
% now plot the solution 
[X,S] = desolve(root);
Nb = alp*L;
figure(2)
plot(X,S(:,2),[Nb,Nb],[0,root],'--','linewidth',2)
 ylabel('c(x) (mM)','fontsize',20)
 xlabel('x (mm)','fontsize',20)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the function transport_eqs(c) solves the initial value problem with
% initial value c and calculates the error for the boundary value condition
% at x=L
function fval = transport_eqs(c)
global D P r c0 alp L N0

[X,S] = desolve(c);
fval = S(end,2)-c0;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,S] = desolve(c)
global D P r c0 alp L N0

%parameters
D=1000;
r=0.05;
c0 = 0.3;
P=0.2;
L=100;
axis([0 L 0 5])
 
alp=0.1;
N0 = 0.3;

xstep = 0.1; % integration step size
x_end = L; % length of  interval
xspan = [0:xstep:x_end];

 % initial data for integration
 y0 = [0;c;0];
   
 stop_cond = odeset('Events',@stopping);   % The stopping conditions for the integration
 %this is needed since for some initial conditions, the integrator cannot
 %integrate all the way to x=L
[X,S] = ode15s(@deRHS,xspan,y0,stop_cond); 
figure(1) 
plot(X,S(:,2))
hold on

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % specify the right hand side of the differential equation system
function s_prime=deRHS(x,s)  % right hand side for ode system
global D P r c0 alp L N0
Nx = N0*(x<=alp*L);
 v=s(1);
 c = s(2);
 cp = s(3);

 vp = 2*P*(c-c0)/r;
 Fc = cp;
 Fcp = (vp*c+v*cp-2*Nx/r)/D;
 
s_prime = [vp,Fc,Fcp]';
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the condition under which to stop the integration 
function [value, isterminal, direction] = stopping(x,y)
value = [y(2);y(2)-5];
isterminal = [1;1];
direction = [0;0];
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the bisection algorithm
function root = bisect(a,b,feval)

ul = a;
fl = feval(ul);
uu = b;
fu = feval(uu);

% we make the assumption, without checking, that 
% fu*fl<0

% if not, the algorithm fails to find a root.
N = 45;  % number of iterates
% the main bisection algorithm
for j = 1:N
uc = (ul+uu)/2;
fc = feval(uc);
ftest = (fc.*fl>0);
ul = ftest.*uc+(1-ftest).*ul;
fl = ftest.*fc + (1-ftest).*fl;

uu = (1-ftest).*uc+ ftest.*uu;
fu = (1-ftest).*fc + ftest.*fu;
end
root = uc;

 
