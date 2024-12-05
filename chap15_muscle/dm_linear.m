% code to reproduce uniform distribution approximation from Math. Biosci. 239 2020:108444
% (Generalized distribution-moment approximation for kinetic theories of muscular contraction)
% GM Donovan
% g.donovan@auckland.ac.nz



function DM_linear
close all
clear all

% for constant stretch (Figs 2 & 3), change sign of v as required
%if(1==0)
%vconst=-10; % +/-
%param.v=@(x)0 + 25*np.sin(50*x);
%end

% for oscillations (Fig 4)
 if(1==1)
 v0=25;
 om=50;
 param.v=@(t) v0*sin(om*t);
 end

% binding functions
f1=43.3;
g1=10;
g2=209;
f=@(x) 0 + (x>0 & x<1).*(f1*x);
g=@(x) g2*(x<=0) + g1*x.*(x>0);
n0=@(x) f(x)./(f(x)+g(x));

% time domain for integration
tmax=0.5; % in s

% parameters into a structure so that they're accessible to the sub-functions
param.h=1;
param.g1=g1;
param.g2=g2;
param.g3=0;
param.f1=f1;

% moments of f. for piecewise linear f you could of course write them explicitly rather than calculate by quadrature. 
xv_orig=linspace(-1,3,1000);
param.b0=trapz(xv_orig,f(xv_orig));
param.b1=trapz(xv_orig,xv_orig.*f(xv_orig));
param.b2=trapz(xv_orig,xv_orig.^2.*f(xv_orig));

cb=colormap('lines').';

% moments of (steady-state) IC
M0_0=trapz(xv_orig,n0(xv_orig));
M1_0=trapz(xv_orig,xv_orig.*n0(xv_orig));
M2_0=trapz(xv_orig,(xv_orig).^2.*n0(xv_orig));

% into a vector
Mbar=[M0_0 M1_0 M2_0];

% integrate
%[tout,yout]=ode45(@(t,y)DMuni_rhs(t,y,param),[0 tvals(end)],Mbar);
[tout,yout]=ode45(@(t,y)DMuni_rhs(t,y,param),[0 tmax],Mbar);

% store
DMuni.t=tout;
DMuni.y=yout;

% plot
for qq=1:3
%figure(qq)
plot(tout,yout(:,qq),'-')
hold on
ylabel(sprintf('M_%d',qq-1))
xlabel('time')
end

% plot the distributions

% set up the subplots
N=10;
sp1=2;
sp2=5;

figure

ttargs=linspace(0,tmax,N);

for qq=1:N
        subplot(sp1,sp2,qq)

        [~,thisind]=min(abs(DMuni.t(:)-ttargs(qq)));
        thisy=DMuni.y(thisind,:);
        M0=thisy(1);
        M1=thisy(2);
        M2=thisy(3);
        a=M1/M0;
        b=(sqrt(3)/M0)*sqrt(M0*M2-M1^2);
        c=(M0^2)/(2*sqrt(3)*sqrt(M0*M2-M1^2));

        plot([a-5*b a-b a-b a+b a+b a+5*b],[0 0 c c 0 0],'color',cb(:,5))
        ax=axis;
        ax(1)=a-4*b;
        ax(2)=a+4*b;
        ax(3)=-0.1;
        ax(4)=1.2;
        axis(ax);
        title(sprintf('t=%0.2f',ttargs(qq)))
end

% labels
subplot(sp1,sp2,1)
ylabel('n')
subplot(sp1,sp2,6)
ylabel('n')
subplot(sp1,sp2,8)
xlabel('x')
subplot(sp1,sp2,10)

end

%% -------------------------------------------------------
function [rhs] = DMuni_rhs(t,y,param)

rhs=0*y; % get the rightshape

M0=y(1);
M1=y(2);
M2=y(3);

% J is implemented in terms of a,b,c instead of M0,M1,M2. (height c between a and b). so get those first
a=M1/M0;
b=(sqrt(3)/M0)*sqrt(M0*M2-M1^2);
c=(M0^2)/(2*sqrt(3)*sqrt(M0*M2-M1^2));

% now the phis
lam=0;
phi0=(param.g2*Juni(lam,a,b,c,-Inf,0) + (param.f1+param.g1)*Juni(lam+1,a,b,c,0,1) + (param.g1+param.g3)*Juni(lam+1,a,b,c,1,Inf) - param.g3*Juni(lam,a,b,c,1,Inf)  );
lam=1;
phi1=(param.g2*Juni(lam,a,b,c,-Inf,0) + (param.f1+param.g1)*Juni(lam+1,a,b,c,0,1) + (param.g1+param.g3)*Juni(lam+1,a,b,c,1,Inf) - param.g3*Juni(lam,a,b,c,1,Inf)  );
lam=2;
phi2=(param.g2*Juni(lam,a,b,c,-Inf,0) + (param.f1+param.g1)*Juni(lam+1,a,b,c,0,1) + (param.g1+param.g3)*Juni(lam+1,a,b,c,1,Inf) - param.g3*Juni(lam,a,b,c,1,Inf)  );

% now the system right-hand side
rhs(1) = param.b0 - phi0;
rhs(2) = param.b1 - phi1 - param.v(t)*y(1);
rhs(3) = param.b2 - phi2 - 2*param.v(t)*y(2);
end

%% -----------------------------------------------------
function z=Juni(lam,a,b,c,alph,bet)
f=@(x) (c/(lam+1))*(x.^(lam+1));
z= f( max( min(a+b,bet),a-b)) - f(min(max(a-b,alph),a+b));
end