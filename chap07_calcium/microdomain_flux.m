% code to solve the microdomain flux equations

%   For Chapter 7, Section 7.6 of
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
global Kb Dbbt Dc De cinf ceinf winf vinf rho 

options = optimset('Display','off');

% parameters
Dc = 1;
De = 1;
Dbbt = 5;
rho = 2;
cinf = 0.5;
ceinf_list = [1 5 10];

n = 1500;
Kb_list = linspace(0.001,10,n);
rho_eff = zeros(n,1);  % preallocation of space, for speed
for j=1:3
    ceinf = ceinf_list(j);
    for i=1:n
        Kb = Kb_list(i);
        winf = getw_v(cinf,Dc);
        vinf = getw_v(ceinf,De);
        c0 = fsolve(@(x)getc0(x),2,options);
        J = rho*(getce0(c0) - c0);
        rho_eff(i) = J/(ceinf - cinf);
    end
    plot(Kb_list,rho_eff)
    hold on
end
xlabel('K_b')
ylabel('\rho_{eff}')
box off
legend('boxoff')
legend('c_{e,\infty} = 1','c_{e,\infty} = 5','c_{e,\infty} = 10')


%% This is the function that must be set to zero, in order to find c0
function out = getc0(c0)
    global Dc De winf vinf
    ce0 = getce0(c0);   % get ce0 as a function of c0. This is done algebraically.

    w0 = getw_v(c0,Dc);
    wp0 = getw_v_prime(c0,Dc);    
    v0 = getw_v(ce0,De);
    vp0 = getw_v_prime(ce0,De);

    out = Dc*(w0 - winf)./(wp0) + De*(v0 - vinf)./(vp0);
end

%%
function ce0 = getce0(c0)
    global Dc winf rho
    w0 = getw_v(c0,Dc);
    wp0 = getw_v_prime(c0,Dc);
    ce0 = c0 + Dc*(w0 - winf)./(wp0*rho);
end

%%
function out = getw_v(x,D)
    global  Dbbt Kb
    out = D.*x + Dbbt.*x./(Kb+x);
end

%%
function out = getw_v_prime(x,D)
    global  Dbbt Kb
    out = D + Dbbt*Kb./((Kb+x).^2);
end



 



