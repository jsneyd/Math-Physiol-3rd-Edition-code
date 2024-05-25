
%  -------------------------------------------------------------------
%
%   3-state model of a L-type calcium channel
%
%   For Chapter 7, Section 7.4.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function L_type_Ca_3state

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

%parameters
p.f=0.85;
p.omega=0.02;
p.g=2;

Vmlist=[-75,  -55, -35, -15, 5,25,45];
for j=1:length(Vmlist)
    Vm=Vmlist(j);

    p.alpha=2*exp(0.012*(Vm-35));
    p.beta = 0.0882*exp(-0.05*(Vm-35));
    c=0;
    p.gamma = 0.44*c;

    p.fbar = p.alpha*p.f/(p.alpha+p.beta);
    p.omegabar = p.omega*(p.beta/2+p.alpha)/(2*p.alpha+p.beta/2);

    dt = 0.01;
    tend=10;
    init = [0.0102,0.9898];
    tspan = [0:dt:tend];
    [t,sol] = ode15s(@(t,x)Ltrhs(t,x,p),tspan,init);
    op=sol(:,1);
    keep(:,j) = op;   % for external plotting

    figure(1)
    plot(t,op )
    hold on
    xlabel('Time (ms)')
    ylabel('Open probability')
end
legend('boxoff')
legend('V=-75','V=-55','V=-35','V=-15','V=5','V=25','V=45')
set(gca,'linewidth',1.5)
box off
hold off

end % of main

%writematrix([tspan', keep],'Ltype.dat')

%%
function out=Ltrhs(t,x,p)

c=(t>3)*2;
p.gamma = 0.44*c;
p.gammabar = p.gamma*(p.beta+2*p.alpha)/(p.alpha+p.beta);
op=x(1);
n3=x(2);
c3=1-n3-op;

opp = -p.g*op+p.fbar*n3;
n3p=-(p.fbar+p.gammabar)*n3 +p.omegabar*c3+p.g*op;

out=[opp;n3p];
end


