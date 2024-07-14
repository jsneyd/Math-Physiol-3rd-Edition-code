
%  -------------------------------------------------------------------
%
%   Phase portraits for leukocyte dynamics. With non-zero alpha.
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
formatSpecF = '%6.2f\n';

%parameters
p.g0=0.2;
p.b0=1;

% choose Icase = 1 or 2;
Icase = 2;

if (Icase ==1)
    alist =[0.5,0.75,1]; % three values of alpha
    xivals=[1.5,2.5,1.7,2.5,1.7];
    bvals = [0.45,.3,2,2,0.5];
    vinit=[0.5,14.5,1.5,1.5,1.0,13.5,2.5];
    tend=[15,30,30,30,30];
end
if(Icase==2)
    xivals=[1.75,2.5,1.6,2.5,2.5];
    bvals = [0.45,.3,2,2,0.5];
    vinit=[0.5,1.4,1.5,1.5,5.1,1.6,5.5];
    tend=[15,40,30,30,300];
    alist=[0.04];% one value of alpha
end

p.a = alist(1);
U0 = p.b0/(1+p.b0);
for j = 1:length(bvals);
    p.b = bvals(j);
    p.xi = xivals(j);
    v=[0:.011:400];
    f = F(p.a*v);
    u1= (p.b0+p.b*v)./((p.b0+p.b*v).*f.*exp(-p.a*v)+ (v.*f+1)) ;
    
    u2=1./(p.xi*f);
    diff = u1-u2;
    ndx =find(diff(1:end-1).*diff(2:end) <=0);
    
    N=1;
    init = [U0 ,vinit(j)];
    if (j==2) N=2;
        init = [U0 ,vinit(j),U0,vinit(6)];
    end
    if (j==5) 
        N=2;
        init = [U0 ,vinit(j),U0,vinit(7)];
    end
    
    %make some phase portraits
    dt=0.01;
    tspan = [0:dt:tend(j)];
    [t,sol] = ode23(@(t,x)rhs(t,x,p),tspan,init);
    
    figure(j)
    plot( sol(:,1) ,sol(:,2),u1,v,'g--',u2,v,'r--','linewidth',2)
    hold on
    plot(U0,0,'k*','linewidth',2)
    if (~isempty(ndx))
        plot(u2(ndx),v(ndx),'k*','linewidth',2)
    end
    
    
    xlabel('U')
    ylabel('V')
    box off
    if(Icase==1)
        if (j==1)
            text(0.36,0.5,'dU/dt=0','fontsize',18)
            text(0.55,1.,'dV/dt=0','fontsize',18)
            annotation('arrow', [.55,.5],[.4,.4])
            annotation('arrow', [.38,.38],[.44,.5])
            plot(0.5,0,'*','linewidth',2)
            title(strcat('A: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
            axis([0.2 0.8 0 3.5])
        end
        if (j==2)
            plot( sol(:,3) ,sol(:,4))
            text(0.34,2,'dU/dt=0','fontsize',18)
            text(0.12,3,'dV/dt=0','fontsize',18)
            annotation('arrow', [.38,.44],[.3,.3])
            annotation('arrow', [.64,.64],[.24,.18])
            axis([.05 .55 0 11])
            title(strcat('B: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
        end
        if (j==3)
            text(0.65,0.85,'dU/dt=0','fontsize',18)
            text(0.47,.55,'dV/dt=0','fontsize',18)
            annotation('arrow', [.54,.49],[.17,.17])
            annotation('arrow', [.25 ,.30],[.4,.4])
            annotation('arrow', [.35,.35],[.14,.2])
            annotation('arrow', [.77,.77],[.57,.51])
            axis([.45 .7 0 2])
            title(strcat('C: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
        end
        if (j==4)
            text(0.63,.66,'dU/dt=0','fontsize',18)
            text(0.35,.75,'dV/dt=0','fontsize',18)
            annotation('arrow', [.19,.24],[.35,.35])
            annotation('arrow', [.82,.82],[.56,.5])
            axis([.3 .7 0 2])
            title(strcat('D: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
        end
        if (j==5)
            text(0.3,4,'dU/dt=0','fontsize',18)
            text(0.25,2.5,'dV/dt=0','fontsize',18)
            annotation('arrow', [.88,.83],[.167,.167])
            annotation('arrow', [.78,.78],[.13,.17])
            axis([.1 .55 0 10])
            title(strcat('E: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
        end
    end
    
    if(Icase==2)
        if (j==1)
            text(0.25,1,'dU/dt=0','fontsize',18)
            text(0.58,1.,'dV/dt=0','fontsize',18)
            annotation('arrow', [.58,.53],[.4,.4])
            annotation('arrow', [.39,.39],[.59,.65])
            plot(0.5,0,'*','linewidth',2)
            title(strcat('A: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
            axis([0 1 0 10])
        end
        if (j==2)
            plot( sol(:,3) ,sol(:,4))
            text(0.3,1.5,'dU/dt=0','fontsize',18)
            text(0.33,0.4,'dV/dt=0','fontsize',18)
            annotation('arrow', [.54,.59],[.3,.3])
            annotation('arrow', [.45,.45],[.75,.82])
            annotation('arrow', [.57,.52],[.7,.7])
            annotation('arrow', [.68,.68],[.27,.21])
            axis([.2 .55 0 2.5])
            title(strcat('B: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
        end
        if (j==3)
            text(0.54,0.6,'dU/dt=0','fontsize',18)
            text(0.625,.6,'dV/dt=0','fontsize',18)
            annotation('arrow', [.63,.58],[.2,.2])
            annotation('arrow', [.52 ,.57],[.7,.7])
            annotation('arrow', [.3,.3],[.18,.24])
            annotation('arrow', [.63,.63],[.86,.8])
            axis([.5 .7 0 2])
            title(strcat('C: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
        end
        if (j==4)
            text(0.55,1.5,'dU/dt=0','fontsize',18)
            text(0.4,1.,'dV/dt=0','fontsize',18)
            annotation('arrow', [.29,.34],[.35,.35])
            annotation('arrow', [.72,.72],[.56,.49])
            axis([.3 .7 0 2])
            title(strcat('D: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
        end
        if (j==5)
            plot( sol(:,3) ,sol(:,4))
            text(0.45,1.,'dU/dt=0','fontsize',18)
            text(0.35,1.,'dV/dt=0','fontsize',18)
            annotation('arrow', [.38,.43],[.25,.25])
            annotation('arrow', [.5,.5],[.27,.2])
            axis([.3 .55 0 8])
            title(strcat('E: \alpha = ',sprintf(formatSpecF,p.a), ', \beta = ',sprintf(formatSpecF,p.b), ', \xi = ',sprintf(formatSpecF,p.xi)),'fontsize',18)
        end
        figure(j)
        hold off
    end
end

% now plot the parameter curves

V=[0.1:.1:1000];

for k = 1:length(alist)
    p.a = alist(k);
    cA=p.a^2*V.^4;
    cB=- V.*(-2*p.a^2*V.^2*p.b0 + p.a*V - 2*exp(p.a*V) + 2);
    cC= p.b0*(p.a^2*V.^2*p.b0 -p.a*V- 1 + exp(p.a*V).*(1 - V.^2*p.a));
    disc=cB.^2-4*cA.*cC;
    b1 = (-cB+sqrt(disc))./(2*cA);
    %b2 = (-cB-sqrt(disc))./(2*cA);
    
    ndx=find(b1>0);
    onebyxi1= p.a*V.*(p.b0 + b1.*V)./(1 - exp(-p.a*V) + p.a*V.^2 + p.a*V.*exp(-p.a*V).*(p.b0 + b1.*V));
    %onebyxi2= a*V.*(p.b0 + b2.*V)./(1 - exp(-a*V) + a*V.^2 + a*V.*exp(-a*V).*(p.b0 + b2.*V));
    
    figure(6)
    ln(k)=plot(1./b1(ndx),1./onebyxi1(ndx));
    hold on
end


xlabel('1/\beta','fontsize',20)
ylabel('\xi','fontsize',20)
text(1,2.5,'D: \{0\}','fontsize',18)
text(0.5,1.3,'C: \{V_p\}','fontsize',18)
if(Icase==2)
    text(1.64,2.12,'\{0,V_p\}','fontsize',16)
end
if(Icase==1)
    text(1.9,1.85,'\{V_p,\infty\}','fontsize',16)
end
text(3,1,'A: \{\infty\}','fontsize',18)
text(3,2.5,'B: \{0,\infty\}','fontsize',18)

figure(6)
plot([0,4],[0,4],'--',[0,4],[1/U0,1/U0],'--' );%,1./bvals,xivals,'r*')  add this if you want to see where the phase portraits are in parameter space
box off
if(Icase==2)
    title(strcat('\alpha = ',sprintf(formatSpecF,p.a)),'fontsize',18)
end
axis([0 4 0 4])

hold off
if (Icase==1)
    legend('boxoff')
    legend([ln(1),ln(2),ln(3)],strcat('\alpha = ',sprintf(formatSpecF,alist(1))),strcat('\alpha = ',sprintf(formatSpecF,alist(2))),strcat('\alpha = ',sprintf(formatSpecF,alist(3))),'fontsize',18,'location','northwest')
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

%%
function f = F(z)
if(z<1.e-4)
    f=1 +  z/2 +  z.^2/12 -  z.^4./720;
else
    f=z./(1-exp(-z));
end



