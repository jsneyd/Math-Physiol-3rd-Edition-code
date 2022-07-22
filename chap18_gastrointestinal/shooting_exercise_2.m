% code to solve the standing gradient water transport equations via
% shooting and bisection

function bisect_function
global L
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

% guess a lower bound for the unknown initial value
a = 0.;

% guess an upper bound for the unknown initial value
b = 2;

%call the function bisect(a,b,@fval)
root = bisect(a,b,@transport_eqs)
% now plot the solution 
[X,S] = desolve(root);
 
figure(2)
plot(X,S(:,1),X,2*sinh(X)./sinh(L),'--','linewidth',2)
 ylabel('u(x) (mM)','fontsize',20)
 xlabel('x (mm)','fontsize',20)
 axis([0 L -2 2])
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the function transport_eqs(c) solves the initial value problem with
% initial value c and calculates the error for the boundary value condition
% at x=L
function fval = transport_eqs(c)


[X,S] = desolve(c);
fval = S(end,1)-2;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,S] = desolve(c)
 global L
L=18;

 

xstep = 0.1; % integration step size
x_end = L; % length of  interval
xspan = [0:xstep:x_end];

 % initial data for integration
 y0 = [0;c];
   
 stop_cond = odeset('Events',@stopping);   % The stopping conditions for the integration
 %this is needed since for some initial conditions, the integrator cannot
 %integrate all the way to x=L
[X,S] = ode15s(@deRHS,xspan,y0,stop_cond); 
figure(1) 
plot(X,S(:,1))
hold on

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % specify the right hand side of the differential equation system
function s_prime=deRHS(x,s)  % right hand side for ode system


 u=s(1);
 up = s(2);
 
 Fu = up;
 Fup = u;
 
s_prime = [Fu Fup]';
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the condition under which to stop the integration 
function [value, isterminal, direction] = stopping(x,y)
value = [y(1)+1;y(1)-10];
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
N = 60;  % number of iterates
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

 

