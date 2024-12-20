
%  -------------------------------------------------------------------
%
%   ODE integrator for the Beeler-Reuter equations.
%
%   For Chapter 12, Section 12.2.3 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  The vector S = [ V m h j d f x Cai] contains                             %
%  V (the membrane potential), m, h, j, d, f, x  gating variables and Cai   %
%  internal Calcium.                                                        %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Beeler_Reuter_ode

close all
clear all
clc

global ENa EK CKi CKe RTbyF t_end gNa0 gK0  gK10 C V0

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);

% set concentrations
RTbyF = 25.8;  % mV
F = 96485.;  % Faraday constant  Coulombs/mol

CNai = 50;      % mM
CNae = 457;     % mM
CKi = 397;      % mM
CKe = 20;       % mM
%CKe = 29; %elevated potassium concentration shows oscillations

% The reversal potentials (mV).
ENa = RTbyF*log(CNae/CNai);  % sodium Nernst potential
EK = RTbyF*log(CKe/CKi);  % potassium Nernst potential
ENa = 50;  % for br


% conductances
gNa0 = 4;
gK0 = 0.8;
gK10 = 0.35;

C = [[0,0,1,-1,-0.1];
    [40,-0.056,0,0,0];
[0.126,-0.25,0,0,0];
[1.7,0,0,1,-0.082];
[0.055,-0.25,0,1,-0.2];
[0.3,0,0,1,-0.1];
[0.095,-0.01,0,1,-0.072];
[0.07,-0.017,0,1,0.05];
[0.012,-0.008,0,1,0.15];
[0.0065,-0.02,0,1,-0.2];
[0.0005,0.0833,0,1,0.057];
[0.0013,-0.06,0,1,-0.04]];


V0=[-47;-72;-77;-22.5;-78;-32;5;-44;-28;-30;-50;-20];

% Steady state initial condition:
s0=[ -84.5738,
	0.011,
	0.9877 ,
	0.9748,
	0.003,
	1,
	0.0056,
	0.0000001782];

t_end = 2000;

tspan = [1: 1:t_end];

[T,S] = ode23s(@BRode,tspan, s0, odeset('maxstep',1));
V = S(:,1);
m = S(:,2);
h = S(:,3);
j = S(:,4);
d = S(:,5);
f = S(:,6);
x = S(:,7);
Ca = S(:,8);

figure(1);
      plot(T,V,'linewidth',2),xlabel('t (ms)','linewidth',2); 
      ylabel('V');
box off
      figure(2)

      plot(T,m,T,h,T,j,'linewidth',2),xlabel('t (ms)','linewidth',2);  
     legend('m','h','j')
     box off

figure(3)
      plot(T,S(:,5),T,S(:,6),T,S(:,7),'linewidth',2),xlabel('t(ms)','linewidth',2); ylabel('d');
      legend('d','f','x')
      box off
     figure(4)
      plot(T,S(:,8),'linewidth',2),xlabel('t (ms)','linewidth',2);
      ylabel('Ca_i');
box off

 for j = 1:length(T)
      currents= IV(S(j,:));
      INa(j) = currents(1);
      IK(j) = currents(2);
      Ix(j) = currents(3);
      Is(j) = currents(4);
 end

figure(5);
     
    plot(T,IK,'r',T,Ix,'b',T,Is,'m','linewidth',2);
    xlabel('t (ms)');
    
    legend('I_K','I_x','I_s')
   
    figure(6)
    plot(T,INa,'g','linewidth',2);
    xlabel('t (ms)');
    ylabel('I_{Na}');

end % of main


%%
function s_prime=BRode(t,s)
    V    = s(1);
    m    = s(2);
    h    = s(3);
    j    = s(4);
    d    = s(5);
    f    = s(6);
    x    = s(7);
    Cai  = s(8);
    gate = [m, h, j, d, f, x]';  % the gating variables


    % the gating variable dynamics

    alpbet=ab(V);
    A=alpbet(1:2:11);
    B=alpbet(2:2:12);

    Inf = A./(A+B);                        %    y_inf = a_y/(a_y + b_y);
    % m = Inf(1);
    Tau = 1./(A+B);
    %    t_y  = 1/(a_y + b_y);

    % Membrane capacitance (micro-Farads/cm^2)
    C_m = 1; %
    % Applied Exterior stimulus

    % only turn on the applied currrent for a short time.
    BCL = 1000;
    if (mod(t,BCL)>10 && mod(t,BCL)<15)
        I = 10;  % stimulus amplitude
    else
        I = 0;
    end

    currents = IV(s);
    INa = currents(1);
    IK = currents(2);
    Ix = currents(3);
    Is = currents(4);
    %
    %
    % The derivatives of V, m, h, j, d, f and x with Cai as well.
    %
    %
    V_prime  = -(INa + IK + Ix + Is - I)/C_m;
    %  gate_prime = (Inf - gate)./Tau;
    gate_prime = A.*(1-gate)-B.*gate;
    Cai_prime = 0.07*(10^(-7) - Cai) - 10^(-7)* Is;
    s_prime = [V_prime ;gate_prime ;Cai_prime];

end

%%
function currents= IV(s)
    global ENa EK gNa0

    V    = s(1);
    m    = s(2);
    h    = s(3);
    j    = s(4);
    d    = s(5);
    f    = s(6);
    x    = s(7);
    Cai  = s(8);

    g_Na = gNa0*m^3*h*j+0.003;                % Na+ conductance

    Ek1=EK;
    % sodium equilibrium potential
    % Currents
    INa = g_Na*(V-ENa);  %sodium current

    %
    % the Ix current
    Ix  = 0.8*x*(exp(0.04*(V+77))-1)/exp(0.04*(V+35));
    % Ix = 0.0084*1.286 *x* (CKi-CKe*exp(-V/25)); %this is the DN current, almost identical

    %the Ik1 current
    IK1= 1.4*(exp(0.04*(V+85))-1)/(exp(0.08*(V+53))+exp(0.04*(V+53))) ...
    +0.07*(V+23)/(1-exp(-0.04*(V+23)));  % the BR current
    IK = IK1  ;

    Esi = -82.3-13.0287*log( Cai);
    Is  = 0.09*d*f*(V-Esi);

    currents = [INa IK Ix Is];

end

%%
function alphabeta=ab(V)
    % there are 5 parameters C for each alpha, beta
    global C V0
    alphabeta = (C(:,1).*exp(C(:,2).*(V-V0))+C(:,3).*(V-V0))./(1+C(:,4).*exp(C(:,5).*(V-V0)));
end

