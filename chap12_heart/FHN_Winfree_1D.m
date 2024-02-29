
%  -------------------------------------------------------------------
%
%  This file looks at solutions of the standard FHN model
%  on a periodic ring with a double stimulus (the Winfree protocol)
%
%   For Chapter 12, Section 12.6.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

clear all
close all
clc

time=cputime;

L = 5; % length of the domain
d = .01; % diffusion coefficient

N=100;
M = 130; % number of time steps
h = L/N;
dt = h/3;


%FHN_parameters;
eps = .05;
alf = .1;
gamma = 2;


  Iamp1 = 0.5;
  Iamp2 = 2;
  t1 = 2;
t2=14;
    scal = d*dt/(2*h^2);

 %FHN_nullclines;
 Vt=0.015*(1:N)-.5;
Vn = -Vt.*(Vt-1).*(Vt-alf);
Wn = Vt/gamma;

   V=zeros(N,1);  % Initial values for V

W = zeros(N,1);  % Initial values for W

    kstep = 20;
Yv = V;
Yw = W;

X = h*(1:N)';
Inx = 1./cosh(5*(X-L/2));
In2x = 1./cosh(4*(X-L/3));
t=0;

% uses Crank Nicolson to solve the diffusion equation

% set up matrix
Ag = scal*ones(N-1,1);

Aexp = -diag([Ag;scal]+[scal;Ag]) +diag(Ag,-1)+diag(Ag,1) ;
Aexp(N,1)=scal;
Aexp(1,N)=scal;

Am = sparse(diag(ones(N,1))-Aexp);

for n=2:M
for k = 1:kstep

 %FHN_dynamics;
 Fw = eps*(V-gamma*W);
Fv  = 20*(-V.*(V-1).*(V-alf)-W);

      W = W + dt*Fw;

Fin = Iamp1/cosh(5*(t-t1))*Inx +Iamp2/cosh(5*(t-t2))*In2x;

V_t =  scal*([V(2:N);V(1)] - 2*V+[V(N);V(1:N-1)]) +dt*(Fv+Fin);

V = Am\(V+V_t);
t = t + dt;
Yv(:,n) = V;
Yw(:,n) = W;

end

figure(1)
  subplot(2,1,1)
  plot(X,V,'r',X,W,'b',X,Fin,'g')
  axis([0 L -.5 1])
 legend('v','w', 'stimulus')
  xlabel('x')
    formatSpecF = '%6.2f\n';

 title(strcat('v, w as functions of x at t=',sprintf(formatSpecF,t)),'fontsize',18)


 w0=2*(alf^2-alf+1)^(3/2)/27;

  subplot(2,1,2)
  plot(V,W,'r',[V(N);V(1)],[W(N);W(1)],'r',Vt,Wn,'k--',Vt,Vn,'b--',Vt,w0,'k--')
  axis([-.5 1 -.1 .2])
   xlabel('v','fontsize',18)
   ylabel('w','fontsize',18)
  title('v-w phase plane','fontsize',[18])

 pause(0.01)
t
end

figure(2)
subplot(2,1,1)
  plot(X,V,'r',X,W,'b')
  axis([0 L -.5 1])
 legend('v','w')
  xlabel('x')
    formatSpecF = '%6.2f\n';

 title('v, w as functions of x','fontsize',18)


 w0=2*(alf^2-alf+1)^(3/2)/27;

  subplot(2,1,2)
  plot(V,W,'r',[V(N);V(1)],[W(N);W(1)],'r',Vt,Wn,'k--',Vt,Vn,'b--',Vt,w0,'k--')
  axis([-.5 1 -.1 .2])
   xlabel('v','fontsize',18)
   ylabel('w','fontsize',18)
  title('v-w phase plane','fontsize',[18])

  %save('onaloop.mat','V','W')
time=cputime-time        % outputs the cpu time taken to solve the equations


