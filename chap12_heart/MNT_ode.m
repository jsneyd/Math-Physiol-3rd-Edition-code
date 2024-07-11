
%  -------------------------------------------------------------------
%
%   ODE integrator for the MNT Purkinje fiber model equations.
%
%   For Chapter 12, Section 12.2.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  The vector S = [ V m h d f q r s x1 x2] contains                             %
%  V (the membrane potential), and gating variables m, h, d, f, q, r, s, x1, x2
                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MNT_ode

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

global gnabar Vna Vcl Vsi C V0 cm Inp
gnabar = 150.;
Vna = 40.;
Vcl=-70;
Vsi=70;

cm = 10.;

%The matrix C has coefficients for alpha beta

C = [[0,0,1,-1,-0.1];
    [40,-0.056,0,0,0];
    [0.0085,-0.184,0 ,0,0];
    [2.5,0,0,1,-0.082];
    [0,0,0.002,-1,-0.1];
    [0.02,-0.089,0,0,0];
[0.000987,-0.04,0,0,0];
[1,0,0,1,-0.087];
[0,0,0.008,-1,-0.1];
[0.08,-0.0888,0,0,0];
[0.00018,-0.04,0,0,0];
[0.02,0,0,1,-0.087];
[0,0,0.001,-1,-0.2];
[5.e-5,-0.067,0,0,0];
[0.0005,0.083,0,1,0.057];
[0.0013,-0.06,0,1,-0.04];
[0.083e-4,0,0,1,-0.2];
[0.0003,-0.06,0,1,-0.04]];

V0=[-47;-72;-71;-10;-40;-40;-60;-26;0;0;-80;-26;-52;-52;-50;-20;-19;-20];

tstep =  0.1;
t_end = 5000;

% specify initial data:
Inp = 0; %.8;
V = -50;
 %initial data
s0 = [V;0.5;0;0.21;0;0;0;1;0.15;0];
%specify the output points
tspan = [0:tstep:t_end];

[T,S] = ode15s(@MNTdeRHS,tspan, s0, odeset('maxstep',1));
V = S(:,1);
m = S(:,2);
h = S(:,3);
d = S(:,4);
f = S(:,5);
q = S(:,6);
r = S(:,7);
s = S(:,8);
x1 = S(:,9);
x2 = S(:,10);

figure(1)
 
plot(T,V,'linewidth',2)
xlabel('t (ms)','fontsize',20)
ylabel('V (mV)','fontsize',20)
box off 
 
figure(2)
plot(T,m,'r',T,h,'b',T,d ,'g','linewidth',2)
legend('m','h','d','fontsize',16)
xlabel('t (ms)','fontsize',20)
box off

figure(3)
plot(T,f ,T,q,T,r,T,s,T,x1,T,x2)
legend('f','q','r','s','x_1','x_2','fontsize',16)
box off
for j = 1:length(T)
    s = S(j,:);
    currents = IV(s);
    INa(j) = currents(1);
    Isi(j) = currents(2);
    IK2(j) = currents(3);
    Ix1(j) = currents(4);
    Ix2(j) = currents(5);
    Icl(j) = currents(6);
    IK1(j) = currents(7);
    Inab(j) = currents(8);
    Iclb(j) = currents(9);

end
figure(4)
plot(T,Isi,T,IK2,T,Ix1,T,Ix2,T,IK1,T,Inab,T,Iclb,'linewidth',2)
legend( 'I_{si}', 'I_{K_2}', 'I_{x_1}', 'I_{x_2}',  'I_{K_1}', 'I_{Na,b}', 'I_{Cl,b}','fontsize',16)
xlabel('t (ms)','fontsize',20)
box off 

figure(5)
plot(T,INa)
legend( 'I_{Na}',  'fontsize',16)
xlabel('t (ms)','fontsize',20)
box off
%

end   % of main


% ---------------------------------------------

%%
%the right hand side for ode simulation:
function s_prime=MNTdeRHS(t,S)
global cm Inp

V = S(1);

alpbet=ab(V);

     currents = IV(S);

      Fv = -(sum(currents) +Inp)/cm;
      Fg = alpbet(1:2:17).*(1.-S(2:10)) - alpbet(2:2:18).*S(2:10);

s_prime = [Fv; Fg];
end

%%
function currents=IV(S)
global gnabar Vna Vsi Vcl
% there are 9 gating variables

V = S(1);
m = S(2);
h = S(3);
d = S(4);
f = S(5);
q = S(6);
r = S(7);
s = S(8);
x1 = S(9);
x2 = S(10);

 gna = gnabar*m.^3*h ;
 INa = gna*(V-Vna);

dpr = 1./(1+exp(-0.15*(V+40)));
Isi= (0.8*d*f+0.04*dpr)*(V-Vsi);

IK2bar = (exp(0.04*(V+110))-1)/(exp(0.08*(V+60))+exp(0.04*(V+60)));
IK2= 2.8*IK2bar*s;

Ix1 = 1.2*x1*(exp(0.04*(V+95))-1)/exp(0.04*(V+45));

Ix2 = x2*(25+0.385*V);

Icl = 2.5*q*r*(V-Vcl);

IK1 = IK2bar+0.2*(V+30)/(1-exp(-0.04*(V+30)));

Inab = 0.105*(V-40);

Iclb = 0.01*(V+70);

currents = [INa, Isi, IK2, Ix1, Ix2, Icl, IK1, Inab, Iclb];
end

%%
function alphabeta=ab(V)
% there are 5 parameters C for each alpha, beta
global C V0
alphabeta = (C(:,1).*exp(C(:,2).*(V-V0))+C(:,3).*(V-V0))./(1+C(:,4).*exp(C(:,5).*(V-V0)));
end
