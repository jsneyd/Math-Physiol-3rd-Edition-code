% code to simulate HH eqn via Method of lines

%  
function Diff_sim
global gnabar gkbar Vna Vk Vl  Veq Tfact N dg sc
 
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 

 % parameters for HH
 Veq = 0; %is the original HH formulation
%Veq = -65; % is the physiological formulation

Tfact = 1; % this is the correction factor for temperatures other that 6.3C
%Tfact = 0.977 % correspnds to 0C

%Tfact = 1.085 % corresponds to 30 C.

gnabar = 120.;
gkbar = 36;
Vna = Tfact*(115+Veq) - Veq;
Vk = Tfact*(-12+Veq)-Veq;
Vl = Tfact*(10.5988+Veq)-Veq;

 % for the continuous cable MOL

 % h = L/N;
 
 L = 40;
 N=401; %number of grid points
 h=L/(N-1);
 dg = 1/(h^2); %Take D = 1  % coupling coefficient
  
 X =[0:N-1]'*h;
 
V = Veq;
m = 0.0529;
h = 0.5961;
n = 0.3177;

 sc = [1;2*ones(N-2,1);1]; % diagonal entries for the diffusion matrix
% initial data
  V0 =  20.1*(1-tanh(X/2)) ;
 
 u0 = [V0;m*ones(N,1);h*ones(N,1);n*ones(N,1)];
 
tstep = 0.05;
t_end = 40;
 
%specify the output points
tspan = [0:tstep:t_end];

[T,S] = ode23(@deRHS,tspan, u0, odeset('maxstep',1));  
 
 
% % Now find the speed:
figure(1)
 mesh(X,T, S(:,1:N))
 xlabel('X')
 ylabel('T')

 thresh = 40;
 % for each X value find the first time the solution crosses the threshold
 for j = 1:N
     jmin = min(find(S(:,j)>=thresh));
     Tc(j) = T(jmin);
 end
 
 xlabel('T')
 ylabel('X')

 p=polyfit(Tc,X,1)
 spest= p(2)+p(1)*Tc;
  speedest = p(1)
  figure(2)
 plot(Tc,X,Tc,spest,'--')
 xlabel('X')
 ylabel('T (ms)')
     formatSpecF = '%5.2f\n';
title(strcat('Speed = ',sprintf(formatSpecF,speedest)),'fontsize',18)
  

 
  
 end
 
%the right hand side for MOL ode simulation:
function s_prime=deRHS(t,s)
global gnabar gkbar Vna Vk Vl  Veq Tfact N dg sc
 % break s up into four parts
 V = s(1:N);
 m = s(N+1:2*N);
 h = s(2*N+1:3*N);
 n = s(3*N+1:4*N);
 
AM=.1*(25.-V)./(exp(.1*(25.-V))-1.);
BM=4.*exp(-V/18.0);
AH=0.07*exp(-V/20.);
BH=1.0./(exp(0.1*(30.-V))+1.);
AN=0.01*(10.001-V)./(exp(.1*(10.001-V))-1.);
BN=0.125*exp(-V/80.);
 
 
Fm=AM.*(1.-m)-BM.*m;
Fh=AH.*(1.-h)-BH.*h;
Fn=AN.*(1.-n)-BN.*n;

currents = IV(s);
INa = currents(:,1);
IK = currents(:,2);
Icl = currents(:,3);

Ion = -INa - IK - Icl;
 
Fv =  dg*(-sc.*V+[0;V(1:end-1)]+[V(2:end);0]) +Ion;
 

s_prime = [Fv ;Fm ;Fh ;Fn];

end

function currents = IV(s)
global gnabar gkbar Vna Vk Vl Veq N
 V = s(1:N);
 m = s(N+1:2*N);
 h = s(2*N+1:3*N); 
 n = s(3*N+1:4*N);
 
gna=gnabar*m.^3.*h;
 
gk=gkbar*n.^4;
Icl=0.3*(V-Vl);
INa = gna.*(V-Vna);
IK = gk.*(V-Vk);

currents = [INa, IK, Icl];

end
