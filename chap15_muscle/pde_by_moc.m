
%  -------------------------------------------------------------------
%
%   Solve the pde u_y +f(u,t,x) u_x = g(u,t,x) by the method of
%   characteristics.
%
%   For Chapter 15, Section 15.12.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------
function pde_by_moc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 
% set the x variable
% spatial grid
dx = 0.01;
x0 = [0:dx:1];
N = length(x0);
%initial data
u =x0.*(1-x0);
x = x0;
t_end = 3;
tstep = 0.25;
%specify the output points
tspan = [0:tstep:t_end];

s0=[x';u'];
figure(1) 
    plot(x,u,'--','linewidth',2)
    xlabel('x','fontsize', 20)
    ylabel('u','fontsize', 20)
    hold on
 
[T,S] = ode23s(@deRHS,tspan, s0 );  

%plot the solution
% each row of S is a parametric representation of u, vs x
Tn = length(T);
N = length(x);
for j = 2:Tn
    figure(1)   
    plot(S(j,1:N),S(j,N+1:2*N),'linewidth',2)
    hold on
end
hold off



%the right hand side for ode simulation:
function s_prime=deRHS(t,s)
     N = length(s);
     x = s(1:N/2);
     u = s(N/2+1:N);
     % specify the MoC equations
     % u_t + Fx u_x = Fu
    Fx = u;
    Fu = zeros(N/2,1);
    s_prime = [Fx; Fu];
