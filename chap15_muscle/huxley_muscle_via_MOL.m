
% to solve u_t + J_x= f(u), with MoL and upwinding (written in conservation
% form, J = v*u), where v is a function of x, t, and u

%    -------------------------------------------------------------------
%
%     Solution of the Huxley model of muscle contraction, using the
%     method of lines. With an oscillatory velocity.
%
%     For Chapter 15, Exercise 15.11 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------


function huxley_muscle_via_MOL
global dx N x f1 g1 g2

%clearvars
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 


par.g1=10;
par.g2=209;
par.f1=43.3;
% put parameters into a structure so that they're accessible to the sub-functions
p.h=1;
 
% spatial grid
dx = 0.01;
x = [-3*p.h:dx:3*p.h];% grid points
N = length(x); % number of points in grid
%initial data:
u0 =f(x,par)./(f(x,par)+g(x,par));

% Use the method of lines (MoL) with upwinding
t_end = 1;  % 
tstep = 0.001;
%specify the output points
tspan = [0:tstep:t_end];
tic
%solve the  MOL differential equations
[T,S] = ode23(@(t,y)deRHS(t,y,par),tspan, u0' );  
% for this problem, ode23 is fastest, ode15s is next and ode23s is slowest
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

end  % of main


%% the right hand side for ode simulation:
function s_prime=deRHS(t,u,par)
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
    Rxn =  (1-u).*f(x',par)-u.*g(x',par);

    s_prime =  Fu + Rxn;
end

%% binding functions
function  out = f(x,par)
    out = 0 + (x>0 & x<1).*(par.f1*x);
end

%%
function  out = g(x,par)
    out = par.g2*(x<=0) + par.g1*x.*(x>0);
end 

%%
function out=v(t)
    % specify the velocity as a function of t
    % in general this can be a function of x, t, and 
    % case 1
    %out = -50;

    % case 2
    out = -50*sin(50*t);
    %out = -20*((t<20).*(t/20) +(t>20).*(2-t/20));
end
