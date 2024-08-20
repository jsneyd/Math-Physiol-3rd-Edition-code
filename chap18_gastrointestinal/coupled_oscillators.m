% 
%    -------------------------------------------------------------------
%
%     Simulate phase equations.
%
%     For Chapter 18, Section 18.5.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function coupled_oscillators
global N delta Omega K x1mG1

close all
clc
set(0,                           ...
    'defaultaxesfontsize', 20,   ...
    'defaultaxeslinewidth', 2.0, ...
    'defaultlinelinewidth', 2.0, ...
    'defaultpatchlinewidth', 0.7);
N = 32; % number of coupled oscillators

x1mG1 = 20;
Omega = -10/31;
num=[1:N];
nat_freq = x1mG1+Omega*(num-1);
nat_freq'
% pick a value of delta
% 
deltalist = [32,18];
for j = 1:length(deltalist)
    delta = deltalist(j);
    
    % coupling matrix
    K = -2*diag(ones(N,1))+diag(ones(N-1,1),-1) + diag(ones(N-1,1),1);
    
    % set up the differential equation solve
    
    tstep = 1; % integration step size
    t_end =2000; % length of  interval
    tspan = [0:tstep:t_end];
    % initial data for integration
    s = zeros(1,N+1);
    [T,S] = ode15s(@deRHS,tspan,s); 
    
    figure(3*(j-1)+1)
        plot(T,S(:,2:end),'linewidth',2)
        xlabel('t','fontsize',20)
        ylabel('\phi_i(t)')
        box off
        formatSpecF = '%6.2f\n';
    figure(3*(j-1)+2)
        plot(S(end,2:N+1)/T(end),'*','linewidth',2)
        xlabel('Oscillator #','fontsize',20)
        ylabel('Average Phase Difference')
        title(strcat('\delta = ',sprintf(formatSpecF,delta)),'fontsize',18)
    box off
    figure(3*(j-1)+3)
        Cs = cumsum(S(end,:))/T(end);
        plot(num,Cs(2:end),'*',num,nat_freq, '--')
        xlabel('Oscillator #','fontsize',20)
        ylabel('Frequency')
        box off
        title(strcat('\delta = ',sprintf(formatSpecF,delta)),'fontsize',18)
end

end % of main

%%
function s_prime=deRHS(t,s)  % right hand side for ode system
    global N delta Omega K x1mG1
    f = sin(s(2:N+1));
    F = Omega + delta*K*f;
    F1 = x1mG1+delta*f(1);
    s_prime = [F1;F];
end
