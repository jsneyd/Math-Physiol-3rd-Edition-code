% This file looks at solutions ofthe standard FHN model
% on a periodic ring following application of a defibrillating stimulus

clear
time=cputime;

L = 5; % length of the domain
d = .01; % diffusion coefficient

N=100;
M = 50; % number of time steps
h = L/N;
dt = h/3;


%FHN_parameters;
eps = .05;
alf = .1;
gamma = 2;

Iampset=[0.85,0.95];
for jamp = 1:2
  Iamp =   Iampset(jamp)
  
  t1 = 2;
t2=14;
    scal = d*dt/(2*h^2);
	
 %FHN_nullclines;
 Vt=0.015*(1:N)-.5;
Vn = -Vt.*(Vt-1).*(Vt-alf);
Wn = Vt/gamma;

 load('onaloop.mat');
     % Initial values for V  and W
 

    kstep = 20;
Yv = V;
Yw = W;

X = h*(1:N)';
 
t=0;

% uses Crank Nicolson to solve the diffusion equation

% set up matrix
Ag = scal*ones(N-1,1);
%Ag(Ng1) = scal*gf1;
%Ag(Ng2) = scal*gf2;


Aexp = -diag([Ag;scal]+[scal;Ag]) +diag(Ag,-1)+diag(Ag,1) ;
Aexp(N,1)=scal;
Aexp(1,N)=scal;

Am = sparse(diag(ones(N,1))-Aexp);

for n=2:M
for k = 1:kstep	

 %FHN_dynamics;
 Fw = eps*(V-gamma*W);
 b= Iamp/cosh(2*(t-t1));
Fv  = 10*(-V.*(V-1).*(V-alf)-W) +b^2*(1+alf-3*V);

      W = W + dt*Fw;
 
V_t =  scal*([V(2:N);V(1)] - 2*V+[V(N);V(1:N-1)]) +dt*(Fv);

V = Am\(V+V_t);
t = t + dt;
Yv(:,n) = V;
Yw(:,n) = W;

end

figure(jamp)
  subplot(2,1,1)
  plot(X,V,'r',X,W,'b' )
  axis([0 L -.5 1])
 legend('v','w' )
  xlabel('x')
    formatSpecF = '%6.2f\n';
 
 title(strcat('v, w as functions of x at t=',sprintf(formatSpecF,t)),'fontsize',18) 
 
 
 w0=2*(alf^2-alf+1)^(3/2)/27;
 
  subplot(2,1,2)

  Vn = -Vt.*(Vt-1).*(Vt-alf) +b^2*(1+alf-3*Vt)/10;;
  plot(V,W,'r',[V(N);V(1)],[W(N);W(1)],'r',Vt,Wn,'k--',Vt,Vn,'b--',Vt,w0,'k--')
  axis([-.5 1 -.1 .2])
   xlabel('v','fontsize',18)
   ylabel('w','fontsize',18)
  title('v-w phase plane','fontsize',[18])
 
 pause(0.01)
t
end

end

   
time=cputime-time        % outputs the cpu time taken to solve the equations


 
