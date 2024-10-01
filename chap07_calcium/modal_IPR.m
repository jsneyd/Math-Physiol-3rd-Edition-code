
%    -------------------------------------------------------------------
%  
%     Compute the steady-state responses and step responses of the simplified
%     modal model.
%  
%     For Chapter 7, Section 7.4.5 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%   
%     Written by James Keener and James Sneyd.
%   
%    ------------------------------------------------------------------- 

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.2, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7); 

global p tau_max  Ktau  Kc Kh  Kp cstim

% Parameter values
tau_max=1000;
Ktau=0.1;
tauP=1;
Kc=0.2;
Kh=0.08;
Kp = 0.2;

keep = [];   % for external plotting

% calculate and plot the steady-state open probability
cc = linspace(0.01,1,2000);
p_list = [0.2,0.4,0.6];
keep = [keep cc'];
for i=1:3
    p = p_list(i);
    h_inf=Kh^4./(Kh^4+cc.^4);
    [~,~,Po]=getstuff(cc,h_inf);
    figure(1)
    semilogx(cc,Po)
    hold on
    keep = [keep Po'];
end
box off
xlim([0,1]); 
xlabel('c (\mu M)');
ylim([0,0.3]); 
ylabel('P_o');
legend('boxoff')
legend('p=0.2 mM','p=0.4 mM','p=0.6mM')


% % Now calculate some responses to step increases in c
init = 1;
p = 0.4;
tspan = linspace(0,30,2000);
keep = [keep tspan'];
cstim_list = [0.4, 0.6, 0.8];
for i = 1:3
    cstim = cstim_list(i);
    [T,Y] = ode15s(@(t,x)oderhs(t,x),tspan,init);
    c = getc(tspan);
    [~,~,Po] = getstuff(c,Y');
    
    % plot solution
    figure(2)
    plot(T,Po)
    box off 
    hold on
    keep = [keep Po'];
end
xlabel('t (s)');
ylabel('P_o');
legend('boxoff')
legend('c_{stim}=0.4 \mu M','c_{stim}=0.6 \mu M','c_{stim}=0.8 \mu M')

%% Define the differential equation here

function out=oderhs(t,h)
    c = getc(t);
    [tauh,h_inf,~] = getstuff(c,h);  
    out = (h_inf-h)/tauh;
end


%%
function [tauh,h_inf,Po]=getstuff(c,h)
    global p tau_max  Ktau  Kc Kh Kp

    phi_c=c.^4./(c.^4+Kc^4);
    phi_p=p^2./(Kp^2+p^2);
    phi_p_down=Kp^2./(Kp^2+p^2);
    h_inf=Kh^4./(Kh^4+c.^4);
    tauh = tau_max*Ktau^4./(Ktau^4+c.^4);
    
    beta = phi_p.*phi_c.*h;
    alpha = phi_p_down.*(1-phi_c.*h_inf);
    Po=beta./(beta + 0.4*(beta+alpha));
    
end
%%
function c=getc(t)
    global cstim
    c = cstim*heaviside(t - 2);
end

