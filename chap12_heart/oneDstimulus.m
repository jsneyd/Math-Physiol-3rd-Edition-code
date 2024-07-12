
%  -------------------------------------------------------------------
%
%   Here we examine the response of a 1D inhomogeneous bidomain cable to a stimulus at the
%   boundaries.
%
%   For Chapter 12, Section 12.6 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function oneDstimulus

clear all
close all
clc

set(0,                           ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0)
global  eps  gam Dscal A  t1 Iamp N Jcol alf

% parameters
eps=0.15;
gam = 0.2;
alf =0.05;

%  %FHN_nullclines;
Vt=[-0.3:0.01:1];
Vn = -10*Vt.*(Vt-1).*(Vt-alf);
Wn = Vt/gam;

N=2000 ; % number of  grid points
K = 10; %number of grid points per cell
Nc = N/K;  %number of cells  (must be an integer)
L=Nc;

dx = 1/K;
x=[1:N]'*dx;
Dscal = 1/dx^2;
ri = 0.1;               % cytoplasmic resistance
re =  0.05;             % extracellular resistance
rglist=[0.1,0.15,0.15]; % gap junctional resistance
Inputlist=[40,45,55];

keep = zeros(N,1);         % for external plotting

for icase = 1:3
    rg = rglist(icase);
    Iamp = Inputlist(icase);
    rc = [ri*ones(1,K-1),rg];
    Rc = [];
    for j = 1:Nc
        Rc = [Rc,rc];
    end
    %remark: we only use N-1 elements of Rc
    R=1./(re+Rc)';

    A=  -spdiags([0;R(1:N-1)],0,N,N) +spdiags(R ,-1,N,N)-spdiags([R(1:N-1);0],0,N,N)  +spdiags([0;R(1:N-1)],1,N,N);
    Jcol = re*([R(1:N-1);0]-[0;R(1:N-1)])/dx;

    t1=2;
    tspan = [0,1,4,6];
    s0 = zeros(2*N,1);
    [T,S] = ode23(@deRHS,tspan, s0);  % Remark ode23s and ode15s are significantly slower for this problem

    for j=1:length(tspan)
        v = S(j,1:N);
        w = S(j,N+1:2*N);

        figure(2*icase-1)
            plot(x,v)
            axis([0 L -0.2 1.2])
            xlabel('Cell Number')
            ylabel('V')
            formatSpecF = '%6.2f\n';
            %
            hold on
            keep = [keep x v'];  % for external plotting

%         figure(2*icase)
%             plot(v,w,Vt,Vn,'--',Vt,Wn,'--')
%             axis([-0.3 1 -1 1.5 ])
    end

    figure(2*icase-1)
        title(strcat('r_g=',sprintf(formatSpecF,rg),', I_{amp} =',sprintf(formatSpecF,Iamp)),'fontsize',18)
        box off
        legend('t=0','t=1','t=4','t=6')
        hold off

end

%writematrix(keep,'oneDstimulus.dat')   % for external plotting

end % of main



%%
%the right hand side for ode simulation:
function s_prime=deRHS(t,u)
    global Dscal A Jcol N gam eps t1 Iamp

    v=u(1:N);
    w=u(N+1:2*N);
    % use method of lines
    b= Iamp*(t<=t1)*(t1^2-t^2);
    vp= Dscal*A*v +  f(v)-w +Jcol*b;
    wp=eps*(v-gam*w);

    s_prime =[vp;wp];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = f(u)
    global alf
    out = 10*u.*(1-u).*(u-alf) ;
end

