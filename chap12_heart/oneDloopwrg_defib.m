% Here we examine the response of a defibrillating stimulus on a one d bidomain  
function onedloop

clear
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)
global  eps  gam Dscal A  t1  t2 Iamp N Jcol alf Iamp1 Iamp2 Inx In2x
tic
% parameters
 eps=0.15;
gam = 0.2;
alf =0.05;
 
 
%  %FHN_nullclines;
 Vt=[-0.3:0.01:1];
Vn = -10*Vt.*(Vt-1).*(Vt-alf);
Wn = Vt/gam;
 

N=2000; % number of  grid points
K = 10; %number of grid points per cell
Nc = N/K;  %number of cells  (must be an integer)
 L=Nc;  

dx = 1/K;
x=[1:N]'*dx;
 
Dscal = 1/dx^2;
ri = 0.1; % cytoplasmic resistance
re =  0.05; % extracellular resistance
 rglist=[0.15,0.15 ]; % gap junctional resistance
Inputlist=[ 610,615];% stimulus amplitude

% read in the initial data
load('onedloopwrg.mat');
     % Initial values for V  and W
 


for icase =1:2
rg=rglist(icase);
Iamp=Inputlist(icase);
Iamp1 = 0;
Iamp2 = 0;

rc = [ri*ones(1,K-1),rg];
Rc = [];
for j = 1:Nc
    Rc = [Rc,rc];
end
 
R=1./(re+Rc)';
  
A=   -spdiags([R(N);R(1:N-1)],0,N,N) +spdiags(R,-1,N,N) -spdiags(R,0,N,N)+spdiags([0;R(1:N-1)],1,N,N);
 A(N,1) = R(N);
 A(1,N) = R(N);
 
Jcol = re*(R-[R(N);R(1:N-1)])/dx;
% 

 Inx = 1./cosh(5*(x-L/2));
In2x = 1./cosh(3*(x-L/3));
t1=2;
t2=36;
 
tspan = [0,1,4,6];
tspan=[0:1:35];
s0 = [v,w];

[T,S] = ode23(@deRHS,tspan, s0);  % Remark ode23s and ode15s are significantly slower for this problem
 
for j=1:length(tspan)
     
 v = S(j,1:N);
 w = S(j,N+1:2*N);
  
 figure(icase)
  formatSpecF = '%6.2f\n';
  subplot(2,1,1)
  % title(strcat( 'I_{amp} =',sprintf(formatSpecF,Iamp)),'fontsize',18) 
 
 plot(x,v)
  axis([0 L -0.3 1.2])
  xlabel('Cell Number')
  ylabel('V')
  
 pause(0.1)
 % 
 
subplot(2,1,2)
plot(v,w,Vt,Vn,'--',Vt,Wn,'--')
 axis([-0.3 1 -1 1.5 ])
 xlabel('v')
 ylabel('w')

end
 
 
 %box off
 % if(icase<3)
 %     legend('boxoff')
 % end
 % 
 %legend('t=0','t=1','t=4','t=6')
 hold off

end
toc

 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the right hand side for method of lines ode simulation:
function s_prime=deRHS(t,u)
global Dscal A Jcol N gam eps t1  t2 Iamp Iamp1 Iamp2 Inx In2x
 
v=u(1:N);
w=u(N+1:2*N);
 b= Iamp/cosh((t-t1)/2);
 
 Fin = Iamp1/cosh(5*(t-t1))*Inx +Iamp2/cosh(5*(t-t2))*In2x;
 vp= Dscal*A*v +  f(v)-w +Fin +Jcol*b;
 wp=eps*(v-gam*w);
 
 s_prime =[vp;wp];
 
 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = f(u)
global alf
 out = 10*u.*(1-u).*(u-alf) ;
end

 