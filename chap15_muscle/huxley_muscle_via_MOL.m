% this code is to solve the Huxley crossbridge mode using upwinding and the method of lines

% to solve u_t + J_x= f(u), with MoL and upwinding (written in conservation
% form, J = v*u), where v is a function of x, t, and u
function pde_MoL_upwind
global dx N x f1 g1 g2

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 


g1=10;
g2=209;
f1=43.3;
% put parameters into a structure so that they're accessible to the sub-functions
p.h=1;
 
% spatial grid
dx = 0.01;
x = [-2*p.h:dx:p.h];% grid points
N = length(x); % number of points in grid
%initial data:
u0 =f(x)./(f(x)+g(x));

% Use the method of lines (MoL) with upwinding
t_end = 0.5;  % 
tstep = 0.005;
%specify the output points
tspan = [0:tstep:t_end];
tic
%solve the  MOL differential equations
[T,S] = ode23(@deRHS,tspan, u0' );  
% for this problem, ode23 is fastest, ode15s is next and ode 23s is slowest
 toc
figure(1)
for j = 1:length(T)
    plot(x,S(j,:),'linewidth',2) 
    
 hold on
end
xlabel('x')
ylabel('h')
hold off

% calculate the force
 F = S*x'*dx;
 figure(2)
 plot(T,F)
xlabel('t')
ylabel('F')
 
figure(3)
 plot(v(T),F)
 xlabel('v')
ylabel('F')

 end
%the right hand side for ode simulation:
function s_prime=deRHS(t,u)
global dx N x
 % by definition, J = v*u, where v = v(x,t,u) in general
% first evaluate v_{j-1/2}
 
xjmh = [x-dx/2,x(N)+dx/2]'; % (x_{j-1/2})
ujmh = ([u;0]+[0;u])/2;  %u_{j-1/2} The zero fill is a ghost point outside the grid

%use xjmh and ujmh to evaluate v_{j-1/2} = v(xjmh,t,ujmh) 

% here is where you evaluate v(x,t,u)
 vjmh = v(t) ;
 
 %this is upwinding:
Jmh = vjmh.*((vjmh>0).*[0;u]+(vjmh<0).*[u;0]);

Fu = (Jmh(1:end-1)-Jmh(2:end))/dx;  %finite difference in x
 
 
s_prime =  (Fu + (1-u).*f(x')-u.*g(x'));
end

% binding functions

 function  out = f(x)
global f1
 
out = 0 + (x>0 & x<1).*(f1*x);
end

function  out = g(x)
global g1 g2

out = g2*(x<=0) + g1*x.*(x>0);
end 

function out=v(t)

% specify the velocity as a function of t
% in general this can be a functioin of x, t, and u

out = -25*sin(50*t);
end

