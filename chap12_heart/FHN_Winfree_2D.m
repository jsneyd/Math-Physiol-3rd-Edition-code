% This file looks at solutions of the standard FHN model
% in a 2d region  with a double stimulus (the Winfree protocol)
function twoD_FHN_Winfree
clear
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)
global alf gam   A eps six siy S1 S2  S3 Iamp1 Iamp2 Iamp3 t1 t2  t3 Nsq

time=cputime;

L = 50; % length/width  of the 2d domain
d = .01; % diffusion coefficient
d=0.05
six=1;
siy = 1;
 
N=90; %number of grid points in each direction 
% must be divisible by 3 and 2
Nsq = N^2;
 
dx = L/N;
 
%FHN_parameters;
eps = .05;
alf = .1;
gam = 2;


X = dx*(1:N)';
Y = X; % a square grid
[X,Y]=meshgrid(dx*(1:N),dx*(1:N));

  Iamp1 = 10;
  Iamp2 = 20;
  Iamp3=20;
X0=L/2;
Y0=L/2;
X1=L/4;
Y1=L/2;
X2=3*L/4;
Y2=L/3;

scale1=0.2;
scale2=1;
scale3 = 0.5;
  %location of stimuli
  Sloc1 = exp(-((X-X0).^2+(Y-Y0).^2)/scale1);
 
  Sloc2 = exp(-((X-X1).^2+(Y-Y1).^2)/scale2);
  Sloc3 = exp(-((X-X2).^2+(Y-Y2).^2)/scale3);
 
  S1=reshape(Sloc1,Nsq,1);
  S2=reshape(Sloc2,Nsq,1);
  S3=reshape(Sloc3,Nsq,1);
 
  t1 = 2;  % time of stimulus 1
  t2=23; %time of stimulus 2
  t3=25; %time of stimulus 3
   
%  %FHN_nullclines;
 Vt=[-0.2:0.01:1];
Vn = -Vt.*(Vt-1).*(Vt-alf);
Wn = Vt/gam;

   V0=zeros(Nsq,1);  % Initial values for V 

W0 = zeros(Nsq,1);  % Initial values for W
 
% uses method of lines to solve the diffusion equation

% set up diffusion matrix

offdiag1 = [ones(N-2,1);2]; % off-diagonal -1
offdiag2 = [2;ones(N-2,1)]; % off-diagonal +1
for j = 1:N-1
    offdiag1 = [offdiag1;0;ones(N-2,1);2];
    offdiag2 = [offdiag2;0;2;ones(N-2,1)];
end
spoffdiag1 = [offdiag1;0];
spoffdiag2 = [0;offdiag2];

offdiag3 = [ones(N*(N-2),1);2*ones(N,1)]; 
spoffdiag3 = offdiag3;
offdiag4 = [2*ones(N,1);ones(N*(N-2),1)];
spoffdiag4 = [zeros(N,1);offdiag4];

A = (-2*(six+siy)*spdiags(ones(Nsq,1),0,Nsq,Nsq)+six*spdiags(spoffdiag2,1,Nsq,Nsq) ...
    + six*spdiags(spoffdiag1,-1,Nsq,Nsq) + siy*spdiags(spoffdiag3,-N,Nsq,Nsq) ...
+siy*spdiags(spoffdiag4,N,Nsq,Nsq))/dx^2;
 
 
 tstep =  1;
t_end = 75;

tspan = [0:tstep:t_end];
 
s0 = [V0;W0 ]

[T,S] = ode23(@deRHS,tspan, s0, odeset('maxstep',1));  
for j = 1:length(T)
    V = reshape(S(j,1:Nsq),N,N);
W=reshape(S(j,Nsq+1:2*Nsq),N,N);
 
figure(1)
 
 %mesh(X,Y,V)
 pcolor(X,Y,V)
 
  xlabel('x')
  ylabel('y')
  zlabel('v')
  axis([0 L 0 L -.25 1])
 
    formatSpecF = '%6.2f\n';
 % 
   title(strcat('t=',sprintf(formatSpecF,T(j))),'fontsize',18) 
 % 

  figure(2)
  plot(V(N/2,:),W(N/2,:),Vt,Vn,'--',Vt,Wn,'--');
  axis([-.25 1,-.2 .2])
 
  pause(0.01)
T(j)
end

time=cputime-time        % outputs the cpu time taken to solve the equations
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the right hand side for ode simulation:
function s_prime=deRHS(t,u)
global Dscal   A gam eps S1 S2 S3 Nsq alf Iamp1 Iamp2 Iamp3 t1 t2 t3
 
V=u(1:Nsq);
 W=u(Nsq+1:2*Nsq);
 %FHN_dynamics;
Fw = eps*(V-gam*W);
Fv  = 10*(-V.*(V-1).*(V-alf)-W);
 
Fin = Iamp1/cosh(5*(t-t1))*S1 +Iamp2/cosh(5*(t-t2))*S2 +Iamp3/cosh(5*(t-t3))*S3;

Vt =   A*V + (Fv+Fin);
 s_prime= [Vt;Fw];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%