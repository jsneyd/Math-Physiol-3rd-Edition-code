%% DDE-Biftool Mackey CN model

clear all
close all
clc

global beta0 phi0 K1 K2


%% Set user-defined functions
beta0      =   8;
phi0      =   0.8;
K2      =   0.36;
K1      =   0.095;
alpha   =   2.4;

gamma_ind=1;
A_ind=2;
tau_ind = [3,4];  % indices of the delays

f=@(N,Ndelay,S,Sdelay,gamma,A,tauN,tauS) [ ...
    -alpha*N + A*phi(Ndelay)*Sdelay; ...
    -phi(N)*S - beta(S)*S + 2*exp(-gamma*tauS)*beta(Sdelay)*Sdelay];
funcs=set_funcs(...
    'sys_rhs',@(xx,p)f(xx(1,1),xx(1,2),xx(2,1),xx(2,3),p(1),p(2),p(3),p(4)),...
    'sys_tau',@()tau_ind);

%% Initial parameters and state. Precomputed, as you can't compute them by hand.
gamma0= 0.0;
A0 = 10.2162;
tauN0=3.5;
tauS0 = 2.8;
N0  = 0.6112;
S0 = 0.4842;

%% Initialization of branch of non-trivial equilibria
contpar=gamma_ind;  % index of the continuation parameter
nontriv_eqs=SetupStst(funcs,'x',[N0;S0],'parameter',[gamma0,A0,tauN0,tauS0],'step',0.01,...
    'contpar',contpar,'max_step',[contpar,0.002],'max_bound',[contpar,1]);

%% Compute and find stability of non-trivial equilibria 
disp('Trivial equilibria');
figure(1);clf
nontriv_eqs=br_contn(funcs,nontriv_eqs,108);
nontriv_eqs=br_stabl(funcs,nontriv_eqs,0,1);
nunst_eqs=GetStability(nontriv_eqs);

%% Compute the periodic branches starting at the Hopf bifurcations
disp('Branch off at first Hopf bifurcation');
ind_hopf=find(nunst_eqs==2,1,'last');
fprintf('Hopf bifurcation near point %d\n',ind_hopf);

[per_orb,suc]=SetupPsol(funcs,nontriv_eqs,ind_hopf,...
    'print_residual_info',0,'intervals',400,'degree',4,...
    'max_bound',[contpar,0.4],'max_step',[contpar,0.05]);
if ~suc
    error('fail',...
        'initialization of periodic orbit failed');
end
figure(1);
hold on
per_orb=br_contn(funcs,per_orb,200);
per_orb=br_stabl(funcs,per_orb,0,1);
nunst_per=GetStability(per_orb,'exclude_trivial',true);

%save('/Users/james/Desktop/Mackey_CN2_results.mat');


%% Read in the results file and plot the max and min of the periodic solutions. Save stuff for external plotting.

%n = size(nunst_eqs(nunst_eqs==2),1); % the number of periodic orbits computed
n = size(per_orb.point,2);
%figure(2)
for i = 1:n
    per(i) = per_orb.point(i).period;
    [min(i),max(i)] = bounds(per_orb.point(i).profile(1,:));
    gam(i) = per_orb.point(i).parameter(1);
end

plot(gam,min,gam,max,'LineWidth',2)

% build and save selected data for graphing in Igor Pro
for i = 1:size(nontriv_eqs.point,2)
    ss(i) = nontriv_eqs.point(i).x(1);
    ss_gam(i) = nontriv_eqs.point(i).parameter(1);
end
writematrix([gam' min' max' per' nunst_per], '/Users/james/Desktop/CN2_per.dat')
writematrix([ss_gam' ss' nunst_eqs], '/Users/james/Desktop/CN2_ss.dat')




%% the phi function
function out=phi(x)
global phi0 K2
    out = phi0*K2./(K2+x);
end
%% the beta function
function out=beta(x)
global beta0 K1
    out = beta0*K1^2./(K1^2+x^2);
end
