% this code is to solve the simplified Huxley-Hai_Murphy smooth muscle crossbridge model 
% using upwinding and the method of lines (MOL)

% to solve u_t + J_x= f(u), with MoL and upwinding (written in conservation
% form, J = v*u), where v is a function of x, t, and u
function pde_MoL_upwind
global   x  

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 

% put parameters into a structure so that they're accessible to the sub-functions

par.f1 = 0.88; 
par.g1 = 0.21;
par.g2 = 4.4;
par.gL1 = 0.01;
par.gL2 = 0.2;
par.h = 1;
par.k1 = 0.35;
par.k2 = 0.1;
par.k5 = 0.35;
par.k6 = 0.1;

par.delx = 1.2*par.h;   % distance between binding sites

 
% spatial grid
par.dx = 0.01;
x = [-4*par.h:par.dx:par.h];% grid points
par.num = length(x) % number of points in spatial grid
%initial data:
nm0 = ones(1,par.num);
nam0 = zeros(1,par.num);
namp0 = zeros(1,par.num);
y0 = [nm0 nam0 namp0];

% Use the method of lines (MoL) with upwinding
t_end = 20;  % 
tstep = 0.05;
%specify the output points
tspan = [0:tstep:t_end];
tic
%solve the  MOL differential equations
[T,S] = ode23(@(t,y)deRHS(t,y,par),tspan,y0 );  
% for this problem, ode23 is fastest, ode15s is next and ode 23s is slowest
 toc
figure(1)
for j = 1:length(T)
    plot(x,S(j,1:par.num),'linewidth',2) 
    
 hold on
end
xlabel('x')
ylabel('N_{M}')
hold off

figure(2)
for j = 1:length(T)
    plot(x,S(j,par.num+1:2*par.num),'linewidth',2) 
    
 hold on
end
xlabel('x')
ylabel('N_{AM}')
hold off

figure(3)
for j = 1:length(T)
    plot(x,S(j,2*par.num+1:3*par.num),'linewidth',2) 
    
 hold on
end
xlabel('x')
ylabel('N_{AM_p}')
hold off


  end
%the right hand side for ode simulation:
function s_prime=deRHS(t,y,par)
global   x
 
nm = y(1:par.num);
nam = y(par.num+1:2*par.num);
namp = y(2*par.num+1:3*par.num);
 nmp = 1 - nm - nam - namp;  % by conservation 


 % by definition, J = v*u, where v = v(x,t,u) in general
 % however, here v=v(t)  so this extra work is not necessary
% first evaluate v_{j-1/2}
 
%xjmh = [x-dx/2,x(N)+dx/2]'; % (x_{j-1/2})
%ujmh = ([u;0]+[0;u])/2;  %u_{j-1/2} The zero fill is a ghost point outside the grid

%use xjmh and ujmh to evaluate v_{j-1/2} = v(xjmh,t,ujmh) 

% here is where you evaluate v(x,t,u)
 vjmh = v(t) ;
 
 %this is the upwinding step:
Jnm = vjmh.*((vjmh>0).*[0;nm]+(vjmh<0).*[nm;0]);
Jnam = vjmh.*((vjmh>0).*[0;nam]+(vjmh<0).*[nam;0]);
Jnamp = vjmh.*((vjmh>0).*[0;namp]+(vjmh<0).*[namp;0]);
 

Fnm = (Jnm(1:end-1)-Jnm(2:end))/par.dx;  %finite difference in x
 Fnam = (Jnam(1:end-1)-Jnam(2:end))/par.dx;  %finite difference in x
Fnamp = (Jnamp(1:end-1)-Jnamp(2:end))/par.dx;  %finite difference in x

% reaction terms
Rxnm =   gL(x',par).*nam -par.k1.*nm + par.k2.*nmp ;
Rxnam = par.k6.*namp-(par.k5+gL(x',par)).*nam;
Rxnamp = par.k5*nam +f(x',par).*nmp -(par.k6+g(x',par)).*namp;
%Rxnmp = par.k1.*nm-(par.k2+f(x,par)).*nmp+g(x,par).*namp;  This for a check
 
s_prime = [Fnm+Rxnm;Fnam+Rxnam;Fnamp+Rxnamp];
 
end

% binding functions

%% ------------
function out = f(x,par)
out = 0 + (x>0 & x<=1).*par.f1.*x;
end

function out = g(x,par)
out = (x<0).*par.g2 + (x>0).*par.g1.*x;
end

function out = gL(x,par)
out = (x<0).*par.gL2 + (x>0).*par.gL1.*x;
end

function out = v(t)
out = -0.6; 
%out = sin(2*t);
end



