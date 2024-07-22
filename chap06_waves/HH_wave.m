%  -------------------------------------------------------------------
%
%   Code to simulate HH eqns via Method of lines
%
%   For Chapter 6, Figure 6.11 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 
%  

function HH_wave

close all
clear all
clc

set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 1.2, ...
'defaultlinelinewidth', 2.0, ...
'defaultpatchlinewidth', 0.7); 

global gnabar gkbar Vna Vk Vl Veq N dg sc gl taum Cm lambda_m


% parameters for HH
Veq = 0;
gnabar = 120.; % mSiemens /cm^2
gkbar = 36;
gl = 0.3;

Cm = 1; % mu F/cm^2
Tfact = 1; % this is the correction factor for temperatures other that 6.3C
%Tfact = 0.977 % correspnds to 0C
%Tfact = 1.085 % corresponds to 30 C.

taum = 1;  % the time constant  in ms
%R=35.4; a=238; taum = 2*R/a   % from Miller/Rinzel paper
lambda_m = 6.5;  % space constant in mm.

Vna = Tfact*(115+Veq) - Veq;
Vk = Tfact*(-12+Veq)-Veq;
Vl = Tfact*(10.5988+Veq)-Veq;

L = 80;
N = 101; %number of grid points
h = L/(N-1);
dg = 1/(h^2); %Take D = 1  % coupling coefficient
X = linspace(0,L,N);

% initial conditions
V = Veq;
m = 0.0529;
h = 0.5961;
n = 0.3177;
V0 = 20.1*(1-tanh(X/2))';
u0 = [V0;m*ones(N,1);h*ones(N,1);n*ones(N,1)];

sc = [1;2*ones(N-2,1);1]; % diagonal entries for the diffusion matrix
tstep = 0.05;
t_end = 60;
tspan = [0:tstep:t_end];

[T,S] = ode15s(@deRHS,tspan, u0, odeset('maxstep',1));  

% % Now find the speed:
figure(1)
    mesh(X,T, S(:,1:N))
    xlabel('X')
    ylabel('T')

% for each X value find the first time the solution crosses the threshold
thresh = 40;
for j = 1:N
    jmin = find(S(:,j)>=thresh,1);
    Tc(j) = T(jmin);
end
% linear fit to find the slope, and thus the wave speed
p=polyfit(Tc,X,1);
spest= p(2)+p(1)*Tc;
speedest = p(1)
figure(2)
    plot(Tc,X,Tc,spest,'--')
    ylabel('X')
    xlabel('T (ms)')
    formatSpecF = '%5.2f\n';
    title(strcat('Speed = ',sprintf(formatSpecF,speedest)),'fontsize',18)

% Plot some wave profiles
figure(3)
    plot(X,S(400,1:N),X,S(600,1:N))

% for external plotting
%keep = [X' S(400,1:N)' S(600,1:N)'];
%writematrix(keep,'hh_wave.dat')

end % of main


%---------------------------------------

%the right hand side for MOL ode simulation:
function s_prime=deRHS(t,s)
    global  N dg sc taum Cm
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
       
    Fm=(AM.*(1.-m)-BM.*m);
    Fh=(AH.*(1.-h)-BH.*h);
    Fn=(AN.*(1.-n)-BN.*n);
    
    [INa,IK,ICl] = IV(V,m,h,n);
    Ion = -INa - IK - ICl;
    
    Fv =  dg*(-sc.*V+[0;V(1:end-1)]+[V(2:end);0])/taum + Ion/Cm ;  
    s_prime = [Fv ;Fm ;Fh ;Fn];

end

%---------------------------------------

function [INa,IK,ICl] = IV(V,m,h,n)
    global gnabar gkbar Vna Vk Vl gl
    gna=gnabar*m.^3.*h;
    gk=gkbar*n.^4;
    ICl=gl*(V-Vl);
    INa = gna.*(V-Vna);
    IK = gk.*(V-Vk);
end
