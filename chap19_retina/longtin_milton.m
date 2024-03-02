
% Code to solve a simple version of the Longtin-Milton model of the pupil
% light reflex. This does not run under Octave, which has not yet implemented
% the delay differential equation solver, dde23.
 
% The Matlab/Octave file used to generate the image in Fig. 25 of Chapter
% 19 
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function longtin_milton
    clear all
    close all
    clc
    
    options = ddeset('Reltol',0.000005);  % To get output at a decent resolution
    tau = 1;   % here is where you set the time delay
    tspan = (0:10);
    sol = dde23(@ddefun, tau, @history, tspan,options);
    area = areafun(sol.y);
    plot(sol.x,sol.y,'b-o',sol.x,area,'r-o')
    xlabel('Time (s)');
    ylabel('Pupil area (mm^2)');
    legend('muscle activity','pupil area')

    % output stuff, for convenience
    % xx = sol.x; yy = sol.y;
    % save('longtin_milton','xx','yy','area')
end

%%
function dydt = ddefun(t,y,Z)
    gamma = 5;
    I = 10;
    phibar = 1;
    taux = 1;
    ylag1 = Z(:,1);
    dydt = (1/taux)*(gamma*capF(log(I*areafun(ylag1)/phibar))-y(1));
end

%%
function s = history(t)
  s = 10*ones(1,1);
end

%% 
function out=capF(x)
out = heaviside(x)*x;
%out = x;
end

%%
function out = areafun(x)
    lam = 30;
    n = 7;   % vary n to get different oscillatory behaviour
    theta = 10;
    out = lam*theta^n./(x.^n + theta^n);
end
