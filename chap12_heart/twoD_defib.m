

%  -------------------------------------------------------------------
%
% This file looks at solutions of the standard FHN model
% in a 2d region following a defibrillating stimulus
%
% All the parameters are set in FHN_Winfree_2D.m, so if you want to change anything (such as the
% resolution) you will need to rerun that code first, otherwise you'll get a mismatch.
% The transfer is done through the file 'doublespiral.mat'.
%
%   For Chapter 12, Sections 12.5.3 and 12.5.6 of
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
'defaultaxeslinewidth', 1.5, ...
'defaultlinelinewidth', 2.0)
global alf gam   A eps six siy  Iamp   t1   Nsq

amplist=[1,1.1];

for ja = 1:length(amplist)
    Iamp=amplist(ja);
    
    load('doublespiral.mat');   % This file loads all the parameters and the initial condition
    S0=[reshape(V,Nsq,1);reshape(W,Nsq,1)];
    X = dx*(1:N)';
    Y = X; % a square grid
    [X,Y]=meshgrid(dx*(1:N),dx*(1:N));
    
    t1 = 2;  % time of stimulus 1
    
    %  %FHN_nullclines;
    Vt=[-0.2:0.01:1];
    Vn = -Vt.*(Vt-1).*(Vt-alf);
    Wn = Vt/gam;
    
    % uses method of lines to solve the diffusion equation
    
    % set up diffusion matrix
    offdiag1 = [ones(N-2,1);2]; % off-diagonal -1
    offdiag2 = [2;ones(N-2,1)]; % off-diagonal +1
    for j = 1:N-1
        offdiag1 = [offdiag1;0;ones(N-2,1);2];
        offdiag2 = [offdiag2;0;2;ones(N-2,1)];
    end
    spoffdiag1 = [offdiag1;0];
    spoffdiag2 = [0;offdiag2];    
    offdiag3 = [ones(N*(N-2),1);2*ones(N,1)]; 
    spoffdiag3 = offdiag3;
    offdiag4 = [2*ones(N,1);ones(N*(N-2),1)];
    spoffdiag4 = [zeros(N,1);offdiag4];
    
    A = (-2*(six+siy)*spdiags(ones(Nsq,1),0,Nsq,Nsq)+six*spdiags(spoffdiag2,1,Nsq,Nsq) ...
    + six*spdiags(spoffdiag1,-1,Nsq,Nsq) + siy*spdiags(spoffdiag3,-N,Nsq,Nsq) ...
    +siy*spdiags(spoffdiag4,N,Nsq,Nsq))/dx^2;
    
    tstep =  0.5;
    t_end = 20;    
    bottom = -0.2;
    top = 1;
    
    tspan = [0:tstep:t_end];
    [T,S] = ode23(@deRHS,tspan, S0, odeset('maxstep',1));  
        
    for j = 1:length(tspan)
        V = reshape(S(j,1:Nsq),N,N);
        W = reshape(S(j,Nsq+1:2*Nsq),N,N);
        
        figure(ja)
            pic = pcolor(X,Y,V);
            pic.LineStyle = 'none';
            clim([bottom top]);      
            xlabel('x')
            ylabel('y')
            zlabel('v')
            axis([0 L 0 L -.3 1])
            formatSpecF = '%6.2f\n';    
            title(strcat('t=',sprintf(formatSpecF,T(j))),'fontsize',18)  

        pause(0.2)  
    end
    
%     keep = zeros(N,N);
%     % save for external plotting
%     jout = [1, 10, 30, 40];
%     for j = 1:length(jout)
%         V = reshape(S(jout(j),1:Nsq),N,N);
%         name = strcat('test',num2str(ja),'_',num2str(j),'.mat');
%         save(name,'V');
%     end
end


%%
%the right hand side for ode simulation:
function s_prime=deRHS(t,u)
    global   A gam eps   Nsq alf Iamp   t1  
    
    V=u(1:Nsq);
    W=u(Nsq+1:2*Nsq);
    %FHN_dynamics;
    Fw = eps*(V-gam*W);
    sc = 4; %  This specifies the length of the stimulus, bigger = shorter
    b= Iamp/cosh( sc*(t-t1));
    Fv  = 10*(-V.*(V-1).*(V-alf)-W) +b^2*(1+alf-3*V);
    
    Vt =   A*V + (Fv );
    s_prime= [Vt;Fw];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%