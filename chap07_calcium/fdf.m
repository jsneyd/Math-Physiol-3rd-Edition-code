
%-------------------------------------------------------------------
%
% The fire-diffuse-fire model.
%
% For Chapter 7, Section 7.8.1 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
% Written by James Keener and James Sneyd
%
%-------------------------------------------------------------------

clear all
close all
clc

set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0);


betalist=[0.5;0.1;0.01]; %;5;10;50;100];
J = length(betalist);
dkappa = 0.001;
eta = [.01:.05:50];
for j = 1:J
    beta = betalist(j);
    for n = 1:1000
        g(n,:) = sqrt(1./(4*pi*n*eta)).*exp(-n./(4*eta) - beta^2*n);
    end    
    G(j,:) = sum(g);
    K(j) = beta;
end

figure(1)
    plot(eta,G,'linewidth',2)
    xlabel('\eta','fontsize',18)
    ylabel('g_\beta(\eta)','fontsize',18)
text(22,0.16,'\beta=0.5','fontsize',18)
text(30,0.54,'\beta=0.1','fontsize',18)
text(22,0.88,'\beta=0.01','fontsize',18)
figure(2)
    plot(G,eta,'linewidth',2)
    ylabel('delay','fontsize',18)
    xlabel('threshold','fontsize',18)
    axis([0 1 0 50])
text(0.15,20,'\beta=0.5','fontsize',18)
text(0.52,20,'\beta=0.1','fontsize',18)
text(0.76,20,'\beta=0.01','fontsize',18)
figure(3)
    semilogy(sqrt(betalist),max(G'),sqrt(betalist),exp(-sqrt(betalist)),'--','linewidth',2)
    xlabel('\beta','fontsize',18)
    ylabel('g_{max}','fontsize',18)
