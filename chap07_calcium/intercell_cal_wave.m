
%  -------------------------------------------------------------------
%
%   Use the method of lines to compute an intercellular calcium wave. This
%   version uses the modal model of the IP3 receptor, but this can be
%   changed. Here, IP3 is also diffusing and has simple reaction terms.
%
%   There are jump conditions specifying
%   the flux across each intercellular boundary.
%
%   For Chapter 7, Section 7.10.2 of
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

p.Vplc = 0.02; % for a regenerating wave
%p.Vplc = 0;     % for no regeneration
p.Kplc = 0.1;
p.kp = 0.1;

p.Dc = 5;       % Ca diffusion coefficient
p.De = 5;       % ER Ca diffusion coefficient
p.DIP = 200;    % IP3 diffusion coefficient
p.FIP = 10;     % intercellular IP3 permeability
p.Fc  = 0;      % intercellular Ca permeability
p.Fe = 0;       % intercellular ER Ca permeability


% Geometry parameters
p.Ncell = 5;        % number of cells
p.N = 500;          % number of spatial grid points. Must be divisible by p.Ncell.
p.nc = p.N/p.Ncell; % number of spatial grid points in each cell.
p.L = 150;          % length of domain
p.delx = p.L/p.N;

p.muc = p.Fc*p.delx/p.Dc;  % used in the cell boundary terms
p.mue = p.Fe*p.delx/p.De;
p.muIP = p.FIP*p.delx/p.DIP;

% set initial data
% all variables start flat except for the IP3 pulse at the left
X = p.delx*(1:p.N)';
c0 = 0.08*ones(p.N,1);
ce0 = 15*ones(p.N,1);
h0 = 0.45*ones(p.N,1);
IP0 = zeros(p.N,1);  
IP0(1:25) = 5;              % IP3 stimulus at left of domain
init = [c0;ce0;h0;IP0];

% specify the output points
t_end = 20;     % total time to run simulation
n_out = 501;    % number of output times
tspan = linspace(0,t_end,n_out);

%%%%%%%%

% integrate the system of odes. It's a DAE because of the jump conditions.
% The zeros of the mass matrix, M, tells you where these jumps occur.

% set up the mass matrix. It needs two zero entries on either side of each
% intecellular boundary, but only for the diffusing species (so not h).
M = eye(4*p.N);
for j = 1:p.Ncell-1
    d = j*p.nc;
    M(d,d) = 0;
    M(d+1,d+1) = 0;
    M(p.N+d,p.N+d) = 0;
    M(p.N+d+1,p.N+d+1) = 0;
    M(3*p.N+d,3*p.N+d) = 0;
    M(3*p.N+d+1,3*p.N+d+1) = 0;
end

options = odeset('Mass',M);
[T,S] = ode15s( @(t,x)pdeRHS(t,x),tspan,init,options); 
c = S(:,1:p.N);             % calcium
ce = S(:,p.N+1:2*p.N);      % ER calcium
IP = S(:,3*p.N+1:4*p.N);    % IP3

% plot a movie  
figure(1)
for i=1:n_out
    plot(X,c(i,:),'linewidth',2)
    ylim([0,1])
    xlim([0,p.L])
    pause(0.1)
end
xlabel('x','fontsize',20)
ylabel('Ca^{++}(x,t)','fontsize',20)

figure(2)
imagesc(X,T,c)
colorbar
xlabel('x (\mu m)')
ylabel('t (s)')





%% the right hand side for pde (MoL) simulation:

function s_prime = pdeRHS(t,s)
global p
 
% There are four variables:  c, ce, h, IP
% only c, ce and IP are diffusing

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

% now put in the intercellular boundaries. A zero term on each side of the
% boundary.
for j=1:p.Ncell-1
    d = j*p.nc;

    Fc(d)   = c(d) - c(d-1) + p.muc*(c(d) - c(d+1));
    Fc(d+1) = c(d+2) - c(d+1) + p.muc*(c(d) - c(d+1));
    
    Fce(d)   = ce(d) - ce(d-1) + p.mue*(ce(d) - ce(d+1));
    Fce(d+1) = ce(d+2) - ce(d+1) + p.mue*(ce(d) - ce(d+1));
    
    FIP(d)   = IP(d) - IP(d-1) + p.muIP*(IP(d) - IP(d+1));
    FIP(d+1) = IP(d+2) - IP(d+1) + p.muIP*(IP(d) - IP(d+1));
end

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
dIPdt = p.Vplc*c.^2./(p.Kplc^2 + c.^2) - p.kp*IP;                          
out = [dcdt;dcedt;dhdt;dIPdt];

end
 
 