
%    -------------------------------------------------------------------
%
%     Solve the Topp model for the pathway to diabetes.
%
%     For Chapter 16, Section 16.7.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------


function topp
close all; clear all; clc;
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0); 

global SI r1 EG0 R0 sigma alpha k d0 r2
SI=1;
EG0=1.44;
R0=864;
sigma=43.2;
alpha=141;
k=432;
d0 = 0.06;
r1 = 1e-3;   
r2 = 0.24e-5;

%% regulated hyperglycemia parameters
 %r1= 1.5;
 %d0 = 0.22;

%% display the steady states and dimensionless parameters

G1 = (r1 + sqrt((r1)^2 - 4*r2*d0))/(2*r2);
G2 = (r1 - sqrt((r1)^2 - 4*r2*d0))/(2*r2);
mu1 = R0/(alpha*EG0);
Ih1 = (R0/G1 - EG0)/SI;
beta1 = (alpha^2 + G1.*G1)./(sigma*G1.*G1) .* k/SI .*(R0./G1 - EG0);
betah1 = k*Ih1*(alpha^2+G1^2)/G1^2;
Ih2 = (R0/G2 - EG0)/SI;
betah2 = k*Ih2*(alpha^2+G2^2)/G2^2;
mu2 = sigma*SI*betah1/(EG0^2);
mu3 = k/EG0;
eps = d0/EG0;
lam1 = r1*alpha/d0;
lam2 = r2*alpha*alpha/d0;

%% fast-slow phase plane
G = linspace(50,1000,1000);
beta = (alpha^2 + G.*G)./(sigma*G.*G) .* k/SI .*(R0./G - EG0);
figure(1)
    plot(beta,G,'--',beta, G1*ones(1,length(beta)),'k--',beta,...
        G2*ones(1,length(beta)),'k--',zeros(1,length(G)),G,'k--', ...
        beta1,G1,'k*',0,R0/EG0,'k*')
    xlim([-5,100])
    xlabel('\beta')
    ylabel('G')
    hold on
%dlmwrite('topp.dat',[beta' G']);

%% solve the ode and plot solution
Idataset = [100,100,10;600,100,12;100,100,17;600,100,17]; %400,100,50;500,100,25;500,100,15;600,100,10];
tspan=linspace(0,100,10000);
for j = 1:4
    [t,sol] = ode15s(@(t,y)rhs(t,y),tspan,Idataset(j,:));
    figure(1)
        plot(sol(:,3),sol(:,1),'LineWidth',2)
        hold on
    figure(2)
        plot(t,sol(:,1))
        hold on
end

figure(1)
    hold off
    figure(2)
    xlabel('t')
    ylabel('G')


end % of main
%% the ODEs
function out = rhs(t,y)
global SI r1 EG0 R0 sigma alpha k d0 r2
G=y(1);
I=y(2);
beta=y(3);
out = [R0 - (EG0+SI*I)*G
       beta*sigma*G*G/(alpha^2+G*G) - k*I
       (-d0+r1*G-r2*G*G)*beta];
end

 