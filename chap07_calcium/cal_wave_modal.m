
%  -------------------------------------------------------------------
%
%   Use the method of lines to compute a traveling calcium wave. This
%   version uses the modal model of the IP3 receptor, but this can be
%   changed. Here, IP3 is also diffusing and has simple reaction terms.
%
%   For Chapter 7, Section 7.6.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 

clear all
close all
clc
set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0);

global p

% parameters

p.tau_max=100;
p.delta=1.5;
p.Ktau=0.1;
p.tauP=1;
p.Kc=0.2;
p.Kh=0.08;
p.kf=10;
 
p.Vserca=0.9;
p.kserca=0.2;
p.Kbar=0.00001957;
p.Kp=0.2;

p.gamma=5.5;
p.Vpm=0.11;
p.Kpm=0.3;
p.alpha0=.0027;
p.Vsocc = .07;
p.Ksocc=8;

p.Vplc = 0.02;
p.kp = 1;

p.Dc = 5;       % Ca diffusion coefficient
p.De = 5;       % ER Ca diffusion coefficient
p.DIP = 200;    % IP3 diffusion coefficient


% integrate the pde
p.N = 600;  % number of spatial grid points
p.L = 60;   % length of domain
p.delx = p.L/p.N;   


% set initial data
X = p.delx*(1:p.N)';
c0 = 0.08*ones(p.N,1);
c0(1:50) = 0.2;  % initial calcium step on left side, to start a wave

% all other variables start flat
ce0 = 15*ones(p.N,1);
h0 = 0.45*ones(p.N,1);
IP0 = (p.Vplc/p.kp)*ones(p.N,1);  % IP3: this is in the range where we expect a solitary pulse
init = [c0;ce0;h0;IP0];

% specify the output points
t_end = 10;     % total time to run simulation
n_out = 101;    % number of output times
tspan = linspace(0,t_end,n_out);

% integrate the system of odes:
[T,S] = ode23( @(t,x)pdeRHS(t,x),tspan,init, odeset('maxstep',1));  

% plot two solutions at different times  
figure(1)
plot(X,S(50,1:p.N),X,S(70,1:p.N),'linewidth',2)
xlabel('x','fontsize',20)
ylabel('Ca^{++}(x,t)','fontsize',20)

%igorpde = [X S(50,1:p.N)' S(70,1:p.N)'];                % for external plotting
%writematrix(igorpde,'igorpde.dat')                      % for external plotting


% Now find the speed
thresh = 0.2;
% for each X value find the first time the solution crosses the threshold
for j = 1:p.N
     jmin = min(find(S(:,j)>=thresh));
     if (jmin) Tc(j) = T(jmin); end   % a check to make sure that jmin exists
end
Xs = X(1:size(Tc,2));
q=polyfit(Tc,Xs,1)
spest= q(2)+q(1)*Tc;
formatSpecF = '%6.2f\n';
figure(4)
plot(Xs,Tc,spest,Tc,'--')
 title(strcat('Velocity = ',sprintf(formatSpecF,q(1)), '\mum s^{-1}'),'fontsize',18)
xlabel('x')
ylabel('t')


%% the right hand side for pde (MoL) simulation:

function s_prime = pdeRHS(t,s)
global p
 
%There are four variables:  c, ce, h, IP
% only IP and C are diffusing

c = s(1:p.N); % calcium
ce = s(p.N+1:2*p.N); % ER calcium
h = s(2*p.N+1:3*p.N); % h
IP = s(3*p.N+1:4*p.N); % IP3

lam_c = p.Dc/p.delx^2;
lam_e = p.De/p.delx^2;
lam_p = p.DIP/p.delx^2;

out=coscrhs(t,s);
sc = [1;2*ones(p.N-2,1);1];

Fc  = lam_c*(-sc.*c + [0;c(1:end-1)]+[c(2:end);0]) + out(1:p.N);
Fce = lam_e*(-sc.*ce + [0;ce(1:end-1)]+[ce(2:end);0])+ out(p.N+1:2*p.N);           
Fh  = out(2*p.N+1:3*p.N);                                                   % reaction only, no diffusion
FIP  = lam_p*(-sc.*IP +[ 0;IP(1:end-1)]+[IP(2:end);0]) + out(3*p.N+1:4*p.N) ;
 
s_prime = [Fc;Fce;Fh;FIP];

end
 

%% the right hand sides of the odes

function out = coscrhs(t,s)
global p

c = s(1:p.N); % calcium
ce = s(p.N+1:2*p.N); % ER calcium
h = s(2*p.N+1:3*p.N); % h
IP = s(3*p.N+1:4*p.N); % IP3

% the equations for the modal model
phi_c=c.^4./(c.^4+p.Kc^4);
phi_p=IP.^2./(p.Kp^2+IP.^2);
phi_p_down=p.Kp^2./(p.Kp^2+IP.^2);

Jpm = p.Vpm*c.^2./(p.Kpm^2+c.^2);
Jin = p.alpha0 + p.Vsocc*(p.Ksocc^4./(p.Ksocc^4+ce.^4));
h_inf=p.Kh^4./(p.Kh^4+c.^4);
tau = p.tau_max*p.Ktau^4./(p.Ktau^4+c.^4);
Jserca=p.Vserca*(c.*c-p.Kbar*ce.*ce)./(c.*c+p.kserca*p.kserca);

beta = phi_p.*phi_c.*h;
alpha = phi_p_down.*(1-phi_c.*h_inf);
Po=beta./(beta + 0.4*(beta+alpha));
Jipr = p.kf*Po.*(ce-c);

dcdt = Jipr-Jserca+p.delta*(Jin-Jpm);
dcedt = p.gamma*(Jserca-Jipr);
dhdt = (h_inf-h)./tau;
dIPdt = p.Vplc - p.kp*IP;                          
out = [dcdt;dcedt;dhdt;dIPdt];

end
 
 