% this code is to solve the Huxley crossbridge mode using  the method of
% characteristics

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
x = [-4*p.h:dx:4*p.h];% grid points
N = length(x); % number of points in grid
%initial data:
u0 =f(x)./(f(x)+g(x));

% Use the method of characteristics
t_end = 1;  % 
tstep = 0.005;
%specify the output points
tspan = [0:tstep:t_end];
tic
%solve the  MOL differential equations
[T,S] = ode23(@deRHS,tspan, [x,u0] );  
% for this problem, ode23 is fastest, ode15s is next and ode 23s is slowest
 toc
 x = S(:,1:N);
 n = S(:,N+1:2*N);

figure(1)
for j = 1:length(T)
   plot(x(j,:),n(j,:),'linewidth',2) 
    
 hold on
end
xlabel('x')
ylabel('h')
hold off

% calculate the force
 F = sum(n.*x*dx,2);
 
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
  x = u(1:N);
  n = u(N+1:2*N);
 
 Fx = -v(t)*ones(N,1);
 
Rxn =  (1-n).*f(x)-n.*g(x);
 size(Fx)
 size(Rxn)
s_prime =  [Fx ; Rxn];
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

out = -20*sin(10*t*pi);
%out = -20*((t<20).*(t/20) +(t>20).*(2-t/20));
end

