
%    -------------------------------------------------------------------
%
%   Model of pulsatile GnRH secretion by Brown et al.
%
%     For Chapter 16, Section 16.2.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function brown
close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
 
par.eps=1/50; par.c=0.2; par.k1=1; par.k2=0.02; par.k3=0.02; par.k4=1; par.rg=2.5;
par.p1=100; par.p2=100; par.p3=0.3; par.rz=1;
formatSpecF = '%6.2f\n';

stimlist = [0.1,0.6,1];
for j = 1:length(stimlist)
    par.stimperiod = stimlist(j);
    par.stimwidth = 0.9*par.stimperiod;
    par.stimheight = 4;
    initial = [0.1 0 0];
    tspan = linspace(0,15,1000);
    [tt,sol]=ode15s(@(t,y)brownfun(t,y,par),tspan,initial);
    figure(2*j-1)
        plot(tt,sol(:,1))
        xlabel('t')
        ylabel('v')
        box off
        title(strcat('Period = ',sprintf(formatSpecF, par.stimperiod)),'fontsize',18)  
        figure(2*j)
        plot(tt,sol(:,3))
        xlabel('t')
        ylabel('z')
        box off
        title(strcat('Period = ',sprintf(formatSpecF, par.stimperiod)),'fontsize',18)  
    %save brown.dat sol -ASCII
end

solve_a_constant(par)
end
%%
function dydt=brownfun(t,y,par)
    v=y(1); g=y(2); z=y(3);
    h = max([v 0]);
    p = par.p1/(1+exp(-par.p2*(h-par.p3)));
    a = getstim(t,par);
    dydt = [(1/par.eps)*(-v*(v-par.c)*(v-1) - par.k1*g + par.k2*a)
            par.k3*a + par.k4*v - par.rg*g
            p - par.rz*z
            ];
end
%%
function a=getstim(t,par)
    a=0;
    if ( rem(t,par.stimperiod)<par.stimwidth) 
        a=par.stimheight;
    end
end

%%
function solve_a_constant(par)
    par.c=0.2; par.k2 = 0.01; par.k3 = 0.05; par.rg = 1;
    % first plot the two sets of nullclines
    a = 0;
    v = linspace(-0.5,1.5,200);
    g1 = v.*(par.c-v).*(v-1) + par.k2*a;
    g2 = par.k3*a + v;
    a = 3;
    g3 = v.*(par.c-v).*(v-1) + par.k2*a;
    g4 = par.k3*a + v;
    figure(10)
    plot(v,g1,'r','HandleVisibility','off')  % don't include in legend
    hold on
    plot(v,g2,'r',v,g3,'b')
    plot(v,g4,'b','HandleVisibility','off')
    ylim([-0.1,0.2])
    hold on
    box off
    
    % now solve for when a is held fixed for various lengths of time, and add
    % to plot. Do this simply by setting a large period, so the stimulus  
    % doesn't repeat
    par.stimperiod = 100;
    par.stimheight = a;
    initial = [0 0 0];
    tspan = linspace(0,10,1000);
    stimlist = [0.1,0.3,1];
    for j = 1:3
        par.stimwidth = stimlist(j);
        [tt,sol]=ode15s(@(t,y)brownfun(t,y,par),tspan,initial);
        plot(sol(:,1),sol(:,2),'LineWidth',2)
    end    
    xlabel('v')
    ylabel('g')
    legend('nullclines for a=0','nullclines for a=3','a=3 for t=0.1','a=3 for t=0.3','a=3 for t=1')
    ax = gca; 
    ax.FontSize = 14;
end
