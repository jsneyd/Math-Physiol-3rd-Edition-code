% ectopic focus oscillations
function ectopy

clear
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)
global alps gam Dscal A x 
% parameters
 
gam = 0.1;
alf =0.15;
b = 0.5;
scale =  1.;
N=50; % number of interior grid points
L = 5; % length of the domain 

dx = L/(N);
x=[1:N]'*dx;
alps=alf +(b-alf)*(exp(-(x/scale).^2));
Dlist = [0;1];
for j = 1: length(Dlist)
    Dl = Dlist(j);
 
Dscal = Dl/dx^2;
A = -2*diag(ones(N,1)) + diag([2;ones(N-2,1)],1) + diag([ ones(N-2,1);2],-1);
 
V0=alps;
tstep =  1;
t_end = 10;

tspan = [0:tstep:t_end];

s0 = [V0 ];

[T,S] = ode23(@deRHS,tspan, s0, odeset('maxstep',1));  
 u = S(end,:);

 
    figure(1)
    plot(x, u )
 hold on
 
figure(2)
plot(x,fp(u))
hold on

% find the eigenvalues
Amat=Dscal*A+diag(fp(u));
max(eig(Amat))
end

% Now the spherical case
A =  diag([-2*ones(N-1,1);-2+2*dx/L]) + diag([ ones(N-1,1)],1) + diag([ ones(N-2,1);2],-1);
 
V0=alps;
tstep = 1;
t_end = 10;

tspan = [0:tstep:t_end];

s0 = [V0 ];

[T,S] = ode23(@deRHSspher,tspan, s0, odeset('maxstep',1));  
 u = S(end,:);
 
    figure(1)
    plot(x, u./x' )
    legend('boxoff')
    legend('uncoupled','1-D symmetric' ,'3-D spherical')
   box off
  xlabel('r')
ylabel('v(r)')
hold off
figure(2)
plot(x,fp(u./x'))
 legend('boxoff')
    legend('uncoupled','1-D symmetric' ,'3-D spherical')
   box off
  xlabel('r')
ylabel('f^{\prime}(v(r))')
hold off

% find the eigenvalues
Amat=Dscal*A+diag(fp(u./x'));
max(eig(Amat))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the right hand side for ode simulation:
function s_prime=deRHS(t,u)
global Dscal   A 
 
 s_prime= Dscal*A*u +  f(u);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the right hand side for ode simulation:
function s_prime=deRHSspher(t,u)
global Dscal   A  x 
 
 s_prime= Dscal*A*u + x.*f(u./x);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = f(u)
global gam alps 
 w=(u-alps)/gam;
out = 10*u.*(1-u).*(u-0.5) -w ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = fp(u)
out = 10*(-3.*u.^2 + 3.*u -0.5);

end
