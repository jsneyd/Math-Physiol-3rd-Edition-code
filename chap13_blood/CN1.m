

%% Be sure to add DDE_Biftool to your Matlab path.

%  -------------------------------------------------------------------
%
%   Use DDE-Biftool to solve the Mackey CN model.
%
%   For Chapter 13, Section 13.1.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

clear all
close all
clc

global delta beta0 th n


% Set user-defined functions
% using gamma and tau as parameters
delta   =   0.08;
beta0   =   1.62;
th      =   3.07;
n       =   4;

gamma_ind=1;
tau_ind=2;

f=@(x,xtau,gamma,tau) 2*exp(-gamma*tau)*beta(xtau)*xtau-(beta(x)+delta)*x;
funcs=set_funcs(...
    'sys_rhs',@(xx,p)f(xx(1,1,:),xx(1,2,:),p(1),p(2)),...
    'sys_tau',@()tau_ind);

% Initial parameters and state
gamma0=0;
tau0=1.39;
x0  = th*(beta0/delta - 1)^(1/n);  % steady state when gamma = 0.

% Initialization of branch of non-trivial equilibria
contpar=gamma_ind;
nontriv_eqs=SetupStst(funcs,'x',x0,'parameter',[gamma0,tau0],'step',0.1,...
    'contpar',contpar,'max_step',[contpar,0.003],'max_bound',[contpar,1]);

% Compute and find stability of non-trivial equilibria 
disp('Trivial equilibria');
figure(1);clf
nontriv_eqs=br_contn(funcs,nontriv_eqs,150); 
nontriv_eqs=br_stabl(funcs,nontriv_eqs,0,1);
nunst_eqs=GetStability(nontriv_eqs);

% Compute the periodic branches starting at the Hopf bifurcations
disp('Branch off at first Hopf bifurcation');
ind_hopf=find(nunst_eqs==2,1,'first');
fprintf('Hopf bifurcation near point %d\n',ind_hopf);

[per_orb,suc]=SetupPsol(funcs,nontriv_eqs,ind_hopf,...
    'print_residual_info',0,'intervals',80,'degree',4,...
    'max_bound',[contpar,80],'max_step',[contpar,0.5]);
if ~suc
    error('fail',...
        'initialization of periodic orbit failed');
end
figure(1);
hold on
per_orb=br_contn(funcs,per_orb,80);
per_orb=br_stabl(funcs,per_orb,0,1);
nunst_per=GetStability(per_orb,'exclude_trivial',true);

save('/Users/james/Desktop/CN1_results.mat');


%% Plot the max and min of the periodic solutions. And other stuff. Save stuff for external plotting.

% DDE-Biftool doesn't plot the orbits this way automatically, so a bit of
% extra work is needed. Not a lot, though.

% Because of saving and reloading the variables, you can run this bit
% independently to save time. Or not.

load('/Users/james/Desktop/CN1_results.mat');

n = size(per_orb.point,2);
for i = 1:n
    per(i) = per_orb.point(i).period;
    [min(i),max(i)] = bounds(per_orb.point(i).profile);
    gam(i) = per_orb.point(i).parameter(1);
end

figure(1)
plot(gam,min,gam,max,'LineWidth',2)
xlabel('\gamma')
ylabel('R')

% plot a typical periodic solution. Save it, if required. Note that the
% horizontal axis is period*mesh, as the mesh is always on [0,1].
xx = per_orb.point(32).mesh*per_orb.point(32).period;
yy = per_orb.point(32).profile;
figure(2)
plot(xx,yy,'LineWidth',2)
xlabel('time')
ylabel('R')
writematrix([xx' yy'],'/Users/james/Desktop/CN_profile.dat')

% % build and save selected data for graphing in Igor Pro
% for i = 1:size(nontriv_eqs.point,2)
%     ss(i) = nontriv_eqs.point(i).x;
%     ss_gam(i) = nontriv_eqs.point(i).parameter(1);
% end
% writematrix([gam' min' max' per'], 'CN_per.dat')
% writematrix([ss_gam' ss' nunst_eqs], 'CN_ss.dat')

%% the beta function
function out=beta(x)
global beta0 th n
    out = beta0*th^n./(th^n+x^n);
end
