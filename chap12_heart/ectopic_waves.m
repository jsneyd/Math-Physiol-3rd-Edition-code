% ectopic focus oscillations
% here we use the metod of lines to compute waves
function ectopic_waves

clear
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)
global alps gam Dscal A x N eps
% parameters

eps = 0.1;
gam = 0.1;
alf =0.;
b = 0.5;
scale = 1.5;
N=200; % number of interior grid points
L = 20; % length of the domain 

dx = L/(N);
x=[1:N]'*dx;
alps=alf +(b-alf)*(exp(-(x/scale).^2));
 
 Dl=1;
Dscal = Dl/dx^2;
% impose a no-flux Neumann boundary condition at both ends
A = -2*diag(ones(N,1)) + diag([2;ones(N-2,1)],1) + diag([ ones(N-2,1);2],-1);
 
V0=zeros(N,1);
W0 = zeros(N,1);
tstep = 0.5;
t_end = 150;

tspan = [0:tstep:t_end];

s0 = [V0;W0 ];

[T,S] = ode23(@deRHS,tspan, s0 );  
[t,X]=meshgrid(T,x);
 
mesh(X,t,S(:,1:N)')
xlabel('x')
ylabel('t')
zlabel('v')


   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the right hand side for MOL ode simulation:
function s_prime=deRHS(t,u)
global Dscal   A N gam  alps eps
 
v=u(1:N);
w=u(N+1:2*N);
 
 vp= Dscal*A*v +  f(v)-w;
 wp=eps*(v-gam*w-alps);
 s_prime =[vp;wp];
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = f(u)
 
out = 10*u.*(1-u).*(u-0.5)  ;
end

 