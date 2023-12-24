%  -------------------------------------------------------------------
%
%   Use the method of lines to compute the traveling wave in the
%   FitzHugh-Nagumo equations. Also see the code fhn_wave.m.
%
%   For Chapter 6, Figure 6.5 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  -------------------------------------------------------------------  

set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 1.2, ...
'defaultlinelinewidth', 2.0, ...
'defaultpatchlinewidth', 0.7); 

% parameters for FHN

p.Iapp=0;
p.alpha=0.2;
p.eps=0.05;
p.gamma=0.5 ;

% for the continuous cable MOL
L = 40;
p.N=500; %number of grid points
h=L/(p.N-1);
p.dg = 1/(h^2); %Take D = 1  % coupling (diffusion) coefficient

X =[0:p.N-1]'*h;

p.sc = [1;2*ones(p.N-2,1);1]; % diagonal entries for the diffusion matrix
% initial data
V0 =   (sech(X));
u0 = [V0;zeros(p.N,1)];
tstep = 0.025;
t_end = 20;

%specify the output points
tspan = [0:tstep:t_end];

[T,S] = ode15s(@(t,x)RHS(t,x,p),tspan,u0);  

%calculate the nullclines
u = [-0.3:.01:1];
w1=(6*u.*(u-p.alpha).*(1-u)+p.Iapp) ;
w2=u/p.gamma;

j=250;

figure(1)
plot(u,w1,'--',u,w2,'--',S(j,1:p.N),S(j,p.N+1:2*p.N),'linewidth',2)
axis([-0.3 1 -0.2 1])
xlabel('v')
ylabel('w')
legend('dv/dt=0','dw/dt=0')


figure(3)
mesh(X,T, S(:,1:p.N))
xlabel('\it x')
ylabel('\it t')
zlabel('\it v')


% Now find the speed:
thresh = 0.5;
% for each X value find the first time the solution crosses the threshold
for j = 1:p.N
     jmin = min(find(S(:,j)>=thresh));
     Tc(j) = T(jmin);
end
q=polyfit(Tc,X,1)
spest= q(2)+q(1)*Tc;

figure(4)
plot(X,Tc,spest,Tc,'--')
xlabel('x')
ylabel('t')
speedest = q(1);


figure(2)
j=250;  % this is an arbitrary value to make a nice graph
plot(X,S(j,1:p.N),X,S(j,p.N+1:2*p.N))
xlabel('x')
legend('v','w')
formatSpecF = '%5.2f\n';
title(strcat('Speed = ',sprintf(formatSpecF,speedest)),'fontsize',18)

% % temporary output for a better plotting program
% vv = S(j,1:p.N);
% ww = S(j,p.N+1:2*p.N);
% XX = X';
% save('igorplot.mat','u','w1','w2','XX','vv','ww')

 
%% the right hand side for MOL ode simulation:
function s_prime=RHS(t,s,p)

% break s up into two parts
V = s(1:p.N);
w = s(p.N+1:2*p.N);

currents = IV(s,p);
FV= currents(:,1);
Fw = currents(:,2);


Fv =  p.dg*(-p.sc.*V+[0;V(1:end-1)]+[V(2:end);0]) +FV;
s_prime = [Fv ;Fw];

end

%%
function currents = IV(s,p)
% calculate the ode dynamics
u = s(1:p.N);
w = s(p.N+1:2*p.N);
Fv = (6*u.*(u-p.alpha).*(1-u)-w +p.Iapp)/p.eps ;
Fw = (u-p.gamma*w) ;

currents = [Fv, Fw];


end
