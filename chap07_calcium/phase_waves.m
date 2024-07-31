%  -------------------------------------------------------------------
%
%   Use the method of lines to compute a traveling calcium wave. This
%   version includes IP3 diffusion, and is used to demonstrate phase waves.
%
%   For Chapter 7, Section 7.9 of
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
    

% Parameter values
par.ct = 2;
par.tau_max=1000;
par.Ktau=0.1;
par.tauP=1;
par.Kc=0.2;
par.Kh=0.08;
par.kf=10;
par.Vserca=0.9;
par.kserca=0.2;
par.Kbar=0.00001957;
par.Kp=0.2;
par.gamma=5.5;


%% integrate the pde

par.N=300;  % number of spatial grid points
par.L=10;
par.h=par.L/par.N;

par.du = 0.3; % diffusion coefficient for IP3
% choose your Ca diffusion coefficient
par.dv =  0.001; 
par.dv =  0.0; 

par.sc = [1;2*ones(par.N-2,1);1];
X = par.h*(1:par.N)';

% choose the initial gradient for P
%P0  = X*(0.1-0.3)/par.L + 0.3;  % initial ip3 concentration has a linear gradient
P0  = 1.5*heaviside(1-X);  % initial ip3 concentration has a blob at one end

C0  = 0.1*ones(par.N,1); % ca concentration
h0 = 0.5*ones(par.N,1);
init = [P0;C0;h0];

% specify the output points
tstep = 0.5; % time between plots
t_end = 150; % total time to run simulation
tspan = [0:tstep:t_end];

% integrate the system of odes:
[T,S] = ode23( @(t,x)pdeRHS(t,x,par),tspan,init, odeset('maxstep',1));  

% % plot the solution at each time step to make a simple movie
% for j = 1:length(T)
%     figure(3)
%     subplot(2,1,1)
%     plot(X,S(j,par.N+1:2*par.N),'linewidth',2)
%     xlabel('x','fontsize',20)
%     ylabel('Ca^{++}(x,t)','fontsize',20)
%     %axis([0 max(X) 0 0.5])
%     subplot(2,1,2)
%     plot(X, S(j,1:par.N),'linewidth',2)
%     xlabel('x','fontsize',20)
%     ylabel('IP_3(x,t)','fontsize',20)
%     %axis([0 max(X) 0 3])
% end
 
figure(5) % plot the final solution
plot(X, S(end,1:par.N),X,S(end,par.N+1:2*par.N),'linewidth',2)
legend('IP_3','Ca^{++}','fontsize',20)
xlabel('x','fontsize',20)

% plot a time course
figure(6)
plot(T,S(:,3*par.N/2),T,S(:,1*par.N/2))
xlabel('time')
ylabel('Ca^{++}(L/2)')

figure(7)
contour(X,T,S(:,par.N+1:2*par.N),'linewidth',2)
xlabel('x')
ylabel('t')
zlabel('Ca^{++}')


%%
%the right hand side for pde (MoL) simulation:
function s_prime=pdeRHS(t,s,par)
 
    %There are four variables:  P, C, h, ce 
    % only P and C are diffusing
    par.h = par.L/par.N;
    scu = par.du/par.h^2;
    scv = par.dv/par.h^2;
     
    P = s(1:par.N);
    C = s(par.N+1:2*par.N);
    %h = s(2*par.N+1:3*par.N);
    %ce = par.gamma*(par.ct - C);
    
    % this is method of lines
    % evaluate the ode part
    out=coscrhs(t,s,par);
     
    Fp=out(1:par.N);
    Fc=out(par.N+1:2*par.N);
    Fh=out(2*par.N+1:3*par.N);
    
    FP = scu*(-par.sc.*P+[0;P(1:end-1)]+[P(2:end);0]) + Fp;
    FC = scv*(-par.sc.*C+[0;C(1:end-1)]+[C(2:end);0]) + Fc ;
     
    s_prime = [FP;FC;Fh];
     
end
   
%%
%the right hand side for ode simulation:

function out=coscrhs(t,s,par)

    P = s(1:par.N); %IP3
    C = s(par.N+1:2*par.N); % calcium
    h = s(2*par.N+1:3*par.N);  % inactivation
    ce = par.gamma*(par.ct - C);
    
    phi_c=C.^4./(C.^4+par.Kc^4);
    phi_p=P.^2./(par.Kp^2+P.^2);
    phi_p_down=par.Kp^2./(par.Kp^2+P.^2);
    h_inf=par.Kh^4./(par.Kh^4+C.^4);
    tauh = par.tau_max*par.Ktau^4./(par.Ktau^4+C.^4);
    
    beta = phi_p.*phi_c.*h;
    alpha = phi_p_down.*(1-phi_c.*h_inf);
    Po=beta./(beta + 0.4*(beta+alpha));
    serca=par.Vserca*(C.*C-par.Kbar*ce.*ce)./(C.*C+par.kserca*par.kserca);
    
    Fp = zeros(par.N,1);  % IP3 is not reacting
    Fc = par.kf*Po.*(ce-C) - serca;
    Fh = (h_inf-h)./tauh; 
    
    out = [Fp;Fc;Fh];
end


 