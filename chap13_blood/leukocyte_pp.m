% 
%  -------------------------------------------------------------------
%
%   Phase portraits for leukocyte dynamics.
%
%   For Chapter 13, Section 13.2.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function leukocyte_dynamics
close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 
global N 

%parameters
p.g0=0.2;
p.b0=1;

%there are four cases:
xivals=[1.6,3.0,1.6,3];
bvals = [0.1,0.1,3,3];
letters=['A','B','C','D'];
vinit=[0.5,2.0,1.5,1.5,2.1];
formatSpecF = '%6.2f\n';
p.a = 0.;
U0 = p.b0/(1+p.b0);

for j = 1:4
    p.b = bvals(j);
    p.xi = xivals(j);
    p.letter=letters(j);
    v=[0:.011:2.5];
    f = F(p.a*v);
    
    %u3=  (p.b0+p.b*v)./((p.b0+p.b*v).*f.*exp(-p.a*u)+ (v*f+1));   
    u1= (p.b0+p.b*v)./((p.b0+p.b*v).*f.*exp(-p.a*v)+ (v.*f+1)) ;    
    u2=1./(p.xi*f);
    diff = u1-u2;
    ndx =find(diff(1:end-1).*diff(2:end) <=0);
    
    N=1;
    init = [U0 ,vinit(j)];
    if (j==2) 
        N=2;
        init = [U0 ,vinit(j),U0,vinit(5)];
    end
    
    
    dt=0.01;
    tend=[5,30,30,30];
    tspan = [0:dt:tend(j)];
    [t,sol] = ode23(@(t,x)rhs(t,x,p),tspan,init);
    
    figure(j)
        plot( sol(:,1) ,sol(:,2),u1,v,'g--',u2,v,'r--','linewidth',2) 
        hold on
        plot(U0,0,'*','linewidth',2)
        if (~isempty(ndx)) 
            plot(u2(ndx),v(ndx),'*','linewidth',2)
        end
        xlabel('U')
        ylabel('V')
        box off
        if (j==1)
            text(0.34,0.5,'dU/dt=0','fontsize',18)
            text(0.55,1.,'dV/dt=0','fontsize',18)
            annotation('arrow', [.78,.73],[.4,.4])
            annotation('arrow', [.2,.2],[.59,.65])
            plot(0.5,0,'*','linewidth',2)
            title(strcat( 'A: \xi = ',sprintf(formatSpecF,p.xi),' \beta = ',sprintf(formatSpecF,p.b)),'fontsize',18)
            axis([.3 .7 0 2])
        end
        if (j==2)
            plot( sol(:,3) ,sol(:,4)) 
            text(0.24,1.8,'dU/dt=0','fontsize',18)
            text(0.34,0.4,'dV/dt=0','fontsize',18)
            annotation('arrow', [.5,.45],[.3,.3])
            annotation('arrow', [.34,.34],[.75,.82])
            annotation('arrow', [.45,.5],[.7,.7])
            annotation('arrow', [.7,.7],[.3,.24])
            axis([.2 .5 0 2.5])
            title(strcat( 'B: \xi = ',sprintf(formatSpecF,p.xi),' \beta = ',sprintf(formatSpecF,p.b)),'fontsize',18)   
        end
        if (j==3)
            text(0.65,0.6,'dU/dt=0','fontsize',18)
            text(0.59,1,'dV/dt=0','fontsize',18)
            annotation('arrow', [.64,.59],[.2,.2])
            annotation('arrow', [.59 ,.64],[.6,.6])
            annotation('arrow', [.4,.4],[.17,.23])
            annotation('arrow', [.87,.87],[.81,.74])
            axis([.5 .7 0 2])
            title(strcat( 'C: \xi = ',sprintf(formatSpecF,p.xi),' \beta = ',sprintf(formatSpecF,p.b)),'fontsize',18)    
        end
        if (j==4)
            text(0.61,1.5,'dU/dt=0','fontsize',18)
            text(0.35,1.,'dV/dt=0','fontsize',18)
            annotation('arrow', [.18,.23],[.35,.35])
            annotation('arrow', [.84,.84],[.56,.5])
            axis([.3 .7 0 2])
            title(strcat( 'D: \xi = ',sprintf(formatSpecF,p.xi),' \beta = ',sprintf(formatSpecF,p.b)),'fontsize',18)
        end
        hold off
end
    
% now plot the parameter curves
% for alpha = 0;

figure(5)
    plot([0,3],[1,4],'--',[0,4],[1/U0,1/U0],'--')
    box off
    xlabel('1/\beta','fontsize',20)
    ylabel('\xi','fontsize',20)
    text(0.7,2.5,'D: \{0\}','fontsize',18)
    text(0.1,1.8,'C: \{V_p\}','fontsize',18)
    text(2.5,1.8,'A: \{\infty\}','fontsize',18)
    text(2.5,2.5,'B: \{0,\infty\}','fontsize',18)

end

%%
function out=rhs(t,x,p)
global N
    out = [];
    for j = 1:N
        u=x(2*j-1); % leukocyte mass
        v=x(2*j); % 
        f = F(p.a*v);
        up = p.g0*((p.b0+p.b*v)*(1-u*f*exp(-p.a*v))-u*(v*f+1));
        vp = v*(1-p.xi*u*f);
        
        out =[out;up;vp];
    end
end

%%
function f = F(z)
    if(z<1.e-4)
        f=1 +  z/2 +  z.^2/12 -  z.^4./720;
    else
        f=z./(1-exp(-z));
    end
end


 