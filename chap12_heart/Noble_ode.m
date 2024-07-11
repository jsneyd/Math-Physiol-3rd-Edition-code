
%  -------------------------------------------------------------------
%
%   ODE integrator for the Noble Purkinje fiber model equations.
%
%   For Chapter 12, Section 12.2.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  The vector S = [ V m h n] contains                             %
%  V (the membrane potential), and gating variables m, h, n.                                                        %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Noble_ode

%this uses a slightly different way of organizing parameters
close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
% The natural period is about 800'
global gnabar Vna gan Van Vk cm C V0 Inp

gnabar = 400.;
Vna = 40.;
gan = 0.14;
Van = 40.;
Vk = -100.;
cm = 12.;

%The matrix C contains parameters for alpha beta's
C=[[0,0,0.1,-1,-1/15];  [0,0,-0.12,-1,0.2];
     [0.17, -0.05 ,0 ,0 ,0];[1,0,0,1,-0.1];
     [0, 0, 0.0001,-1,-0.1];[0.002,-0.0125,0,0,0]];

V0=[-48;-8;-90;-42;-50;-90];  %Nernst potentials of the currents

tstep = .2;
t_end = 3500;

% specify initial data:
Inp = 0; %.8;
V = -80;
m=0.1;
h = 0.8;
n = 0.3;
%specify the output points
tspan = [0:tstep:t_end];

s0 = [V;m;h;n];

[T,S] = ode23(@NobledeRHS,tspan, s0, odeset('maxstep',1));
V = S(:,1);
m = S(:,2);
h = S(:,3);
n = S(:,4);

figure(1)
 
plot(T,V,'linewidth',2)
xlabel('t (ms)','fontsize',20)
ylabel('V (mV)','fontsize',20)
box off
figure(2)
plot(T,m,'r',T,h,'b',T,n,'g','linewidth',2)
legend('m','h','n','fontsize',16)
xlabel('t (ms)','fontsize',20)
box off
for j = 1:length(T)
    s = S(j,:);
    currents = IV(s);
    INa(j) = currents(1);
    IK(j) = currents(2);
    Icl(j) = currents(3);
end
figure(3)
plot(T,INa,T,IK,T,Icl,'linewidth',2)
legend('I_{Na}','I_K','I_l','fontsize',16)
xlabel('t (ms)','fontsize',20)
box off
end % of main

% ------------------------------------------


%%
%the right hand side for ode simulation:
function s_prime=NobledeRHS(t,s)
global gnabar Vna gan Van Vk C Inp cm

V = s(1);
m = s(2);
n = s(4);
h = s(3);

alphabeta=ab(V);
currents = IV(s);
INa = currents(1);
IK = currents(2);
Ian = currents(3);

Fv = -(INa+IK+Ian +Inp)/cm;
Fg = alphabeta(1:2:5).*(1-s(2:4))-alphabeta(2:2:6).*s(2:4);

s_prime = [Fv ; Fg];
end


%%
function currents=IV(s)
global gnabar Vna gan Van Vk cm Inp

V = s(1);
m = s(2);
n = s(4);
h = s(3);

gna = gnabar*m.^3*h ;
gk1 = 1.2*exp(-(V+90.)/50.) + 0.015*exp((V+90.)/60.);
gk2 = 1.2*n.^4;
INa = gna*(V-Vna);
Ik = (gk1+gk2)*(V-Vk);
Ian = gan*(V-Van);

currents = [INa, Ik, Ian];
end


%%
function alphabeta=ab(V)
% there are 5 parameters C for each alpha, beta
global C V0
alphabeta = (C(:,1).*exp(C(:,2).*(V-V0))+C(:,3).*(V-V0))./(1+C(:,4).*exp(C(:,5).*(V-V0)));
end
