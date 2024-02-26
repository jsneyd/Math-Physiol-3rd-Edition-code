
%  -------------------------------------------------------------------
%
%   Plot circle maps F(T_{n+1} ) = G(T_n).
%
%   For Chapter 12, Section 12.5.2 of
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
   'defaultlinelinewidth', 2.0)

global gam del

% set parameters 
del = 1;
gamlist=[0.75,0.694,0.67,0.55];  % a list of parameter values to use

for ngam = 1:length(gamlist)
    clear t1 t2 t3 tk Gkj Tk
    gam=gamlist(ngam);
    
    delt=0.001; % This determines the resolution 
    t=[0:delt:6];  
    F = evalF(t);    
    G=evalG(t);
    
    % iterate the map
    
    tk(1) = 0.;
    N=0;
    j=1;
    while N==0
        tj = tk(j);
        Gkj=evalG(tj);
        Gk(j) = Gkj;
        if(isempty(find(F>=Gkj)))
            N=1;
        else
            j=j+1;
            ndx = min(find(F>=Gkj));
            tk(j) =t(ndx);
        end
    end 
    N=j;
    gk(1) = 0;
    for j = 1:N
        gk(2*j) = Gk(j);
        gk(2*j+1) = Gk(j);
        Tk(2*j-1) = tk(j);
        Tk(2*j) = tk(j);
    end
    
    %this shows an example of the map
    figure(2*ngam-1)
        plot(t,G,t,F,Tk,gk(1:2*N),'k--')
        axis([0 4 0 20])
        xlabel('t','fontsize',20)
        text(3.65,8,'F(t)','fontsize',18)
        text(3.65,16.5,'G(t)','fontsize',18)
        text(0.3,3.75,'t_1','fontsize',18)
        text(3.1,12.3,'t_3','fontsize',18)
        text(1.3,5.95,'t_2','fontsize',18)
        formatSpecF = '%6.2f\n';
        title(strcat('\gammaT = ',sprintf(formatSpecF,gam)),'fontsize',18)
    
    % plot some circle maps (on the circle)    
    t=[0:delt:1];
    F = evalF(t); 
    [M,I]=max(F);
    ti=[0:delt:t(I)];
    Gkj=evalG(ti); 
    t1=[0:delt:3];  %this must be large enough to get the whole range, but not too large
    F1=evalF(t1);
    for j=1:I
        j1=min(find(F1>=Gkj(j)));  
        t2(j)= t1(j1);
    end
    t3=mod(t2,1);
        
    % track a trajectory
    tk(1)=0.35; % starting value
    t=[0:delt:15];
    F = evalF(t); 
    N=0;
    j=1;
    while N==0
        tj = tk(j);
        Gkj=evalG(tj);
        if(isempty(find(F>=Gkj)))
            N=1;
        else
            j=j+1;
            ndx = min(find(F>=Gkj));
            tk(j) =t(ndx);
        end
    end
    N=j;
    
    for j = 1:N        
        Tk(2*j-1) = tk(j);
        Tk(2*j) = tk(j);
    end
    Istop=max(find(t3(2:end)<t3(1:end-1))); % find the discontinuity/truncation
    figure(2*ngam)
    if(isempty(Istop))
        plot(ti,t3,'r',t3,t3,'b--',mod(Tk(1:end-1),1),mod(Tk(2:end),1),'k--')
        axis([min(t3) max(t3) min(t3) max(t3)])
    else
        plot(ti(1:Istop),t3(1:Istop),'r',ti(Istop+1:end),t3(Istop+1:end),'r',t3,t3,'b--',mod(Tk(1:end-1),1),mod(Tk(2:end),1),'k--')
        axis([min(t3) max(t3) min(t3) max(t3)]) 
    end
    
    xlabel('\Psi_n')
    ylabel('\Psi_{n+1}')
    box off
    formatSpecF = '%6.2f\n';
    title(strcat('\gammaT = ',sprintf(formatSpecF,gam)),'fontsize',18)
end

%% here is where the functions F and G are specified
function F = evalF(t)
    global gam
    F = sin(pi*t).^4.*exp(gam*t);
end

function G = evalG(t)
    global gam del 
    G = evalF(t)+del*exp(gam*t);
end


 