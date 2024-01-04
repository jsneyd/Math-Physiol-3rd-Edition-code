% to solve u_t + J_x= 0, with MoL and upwinding (written in conservation
% form, J = v*u
function pde_MoL_upwind
global dx N x

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 

% spatial grid
dx = 0.01;
x = [-3:dx:3];% grid points
N = length(x); % number of points in grid
%initial data:
u =exp(-x.^2);
u = zeros(1,N);

% try method of lines (MoL) with upwinding
t_end = 10;  % 
tstep = 0.1;
%specify the output points
tspan = [0:tstep:t_end];

u0 = u'; % must be a column vector
%solve the  MOL differential equations
[T,S] = ode23s(@deRHS,tspan, u0 );  

 
figure(1)
for j = 1:length(T)
  
    %plot(x,S(j,:)); %, x, x.*exp(T(j)).*(1- x*exp(T(j))),'--','linewidth',2)
 mesh(S)
 %legend('numerical','exact','fontsize',16)
 hold on
end
hold off
 
xlabel('x','fontsize',20)
ylabel('u','fontsize',20)

%the right hand side for ode simulation:
function s_prime=deRHS(t,u)
global dx N x

% by definition, J = v*u, where v = v(x,t,u) in general
% first evaluate v_{j-1/2}
 
xjmh = [x-dx/2,x(N)+dx/2]'; % (x_{j-1/2})
ujmh = ([u;0]+[0;u])/2;  %u_{j-1/2} The zero fill is a ghost point outside the grid

%use xjmh and ujmh to evaluate v_{j-1/2} = v(xjmh,t,ujmh) 

% for this example, v = -x
%vjmh = -  xjmh; 
vjmh=-2*sin(3*t);
 %this is upwinding: 
Jmh = vjmh.*((vjmh>0).*[0;u]+(vjmh<0).*[u;0]);

Fu = (Jmh(1:end-1)-Jmh(2:end))/dx;  %finite difference in x
 
 rxn = (x>0&x<1)'.*(1-u)- 0.5*u;
s_prime =  Fu+rxn;