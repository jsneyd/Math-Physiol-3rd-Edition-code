% using piecewise linear FHN
% here we use the method of lines to compute waves
function ectopic_waves

clear
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)
global alf gam Dscal A x N eps dx Ic 
% parameters

eps = 0.1;
gam = 0.1;
 
alf =0.15;
Ic = 150; % current input  
 
N=100; % number of interior grid points
L = 10; % length of the domain 

dx = L/(N);
x=[1:N]'*dx;
 
 Dl=1;  %diffusion coefficient
Dscal = Dl/dx^2;
% impose a  Neumann boundary condition at both ends
A = -2*diag(ones(N,1)) + diag([2;ones(N-2,1)],1) + diag([ ones(N-2,1);2],-1);
 %initial data
V0=ones(N,1)*alf/(1+gam);;
W0 = zeros(N,1);
tstep = 0.5;
t_end = 200;

tspan = [0:tstep:t_end];

s0 = [V0;W0 ];

[T,S] = ode23(@deRHS,tspan, s0 );  
[t,X]=meshgrid(T,x);
 figure(1)
mesh(X,t,S(:,1:N)')
xlabel('x')
ylabel('t')
zlabel('v')
 
figure(2)
plot(x,S(end,1:N),x,0.25*ones(1,N),'--')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the right hand side for MOL ode simulation:
function s_prime=deRHS(t,u)
global Dscal   A N gam  alf eps Ic dx
 
v=u(1:N);
w=u(N+1:2*N);
 
 vp= Dscal*A*v +  f(v)-w;
 vp(1) = vp(1) +2*Ic/dx;
 wp=eps*(v-gam*w-alf);
 s_prime =[vp;wp];
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = f(u)
 
out = -u.*(u<=1/4) + (u>1/4&u<3/4).*(u-1/2) + (u>=3/4).*(1-u);
end

 