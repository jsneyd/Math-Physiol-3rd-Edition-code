%
%  -------------------------------------------------------------------
%
%   Simulation of the Tyson-Novak model of the fission yeast cell cycle.
%
%   For Chapter 10, Section 10.4.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function fission_yeast

clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);

global K
global k1 k2p k2pp k2ppp
global k3p k3pp k4p  k4
global k5p k5pp k6
global J3 J4 J5
global k7 k8 J7 J8
global k9 k10 J9 J10;
global k11 k12 k12p k12pp Kdiss
global k13p k13pp k14
global k15p k15pp k16p k16pp J15 J16
global mu
global kweep   kweepp   Vawee   Viweep Viweepp Jawee Jiwee
global k25p k25pp Va25p Va25pp Vi25 Ja25 Ji25

formatSpecF = '%6.2f\n';
K=100;
k1=.03; k2p=.03; k2pp=1.; k2ppp=.1;
k3p=1.; k3pp=10.; k4p=2.; k4=35.;
k5p=.005; k5pp=.3; k6=.1;
J3=.01; J4=.01; J5=.3;
k7=1.; k8=.25; J7=.001; J8=.001;
k9=.1; k10=.04;  J9=0.01;  J10=0.01;
k11=.1; k12=0.01; k12p=1.; k12pp=3.; Kdiss=.001;
k13p=0.; k13pp=.1; k14=.1;
k15p=1.5; k15pp=0.; k16p=1.; k16pp=2.; J15=.01; J16=.01;
mu=.005;
kweep=0.15;
kweepp=1.3;

 % to simulate the wee1 mutant set kweepp=0.3;
 kweepp=0.3;
Vawee=0.25;  Viweep=0.;  Viweepp=1.;  Jawee=0.01;  Jiwee=0.01;
k25p=0.05;  k25pp=5;  Va25p=0.;  Va25pp=1.;  Vi25=0.25;  Ja25=0.01;  Ji25=0.01;
kweepplist=[1.3,0.3];
tlim =[0,200];
ylim=[2,1];

for icase=1:2
kweepp=kweepplist(icase);

npts=10; delt=0.2;
numtimesteps=2500; test=0;

y0(1)=0.0029127;         % CycBT
y0(2)=0.001591 ;         % pMPF
y0(3)=0.99997;           % Ste9
y0(4)=0.05 ;             % Slp1T
y0(5)=0.0 ;              % Slp1
y0(6)=0. ;               % IEP
y0(7)=8.5145;            % Rum1T
y0(8)=0.0017441;         % SK
y0(9)=1;                 % Mass
y0(10)=0;                % x (MPF)

for loop=1:numtimesteps
    keep(loop,1)=(loop-1)*delt;         % Keep the times here
    keep(loop,2:11)=y0;                 % Keep the solutions here
    tspan=linspace(0,delt,npts);
    [T,Y] = ode23(@cilibertoodes,tspan,y0);
    if (Y(npts,10)>0.4) test=1; end
    if ((Y(npts,10)<0.05) & (test==1))  % Here is the threshold for x (which tracks MPF)
        Y(npts,9)=Y(npts,9)/2.0;        % Reset the mass for cell division
        test=0;
    end
    y0=Y(npts,1:10);
end

% y0(9)=1.9;
% tspan=linspace(0,400,400);
% [T,Y]=ode23tb(@cilibertoodes,tspan,y0);
% plot(T,Y(:,10))
figure( icase)

 plot(keep(:,1),keep(:,2),keep(:,1),keep(:,10),keep(:,1),keep(:,4),keep(:,1),keep(:,11))
        % Plot the mass,  MPF (or dum, the tracking variable),CDC_T and Ste9
xlabel('t (min)')
box off
 

legend('[Cdc13_T','m','[Ste9]','[MPF]')
 title(strcat('Kwee1 = ',sprintf(formatSpecF,kweepp)),'fontsize',18)
 axis([tlim(icase) 500 0 2])


end



%save cili.dat keep -ascii

end % of main



%%

function dy=cilibertoodes(t,y)

global K
global k1 k2p k2pp k2ppp
global k3p k3pp k4p  k4
global k5p k5pp k6
global J3 J4 J5
global k7 k8 J7 J8
global k9 k10 J9 J10;
global k11 k12 k12p k12pp Kdiss
global k13p k13pp k14
global k15p k15pp k16p k16pp J15 J16
global mu
global kweep   kweepp   Vawee   Viweep Viweepp Jawee Jiwee
global k25p k25pp Va25p Va25pp Vi25 Ja25 Ji25

dy=zeros(10,1);

Cdc13T=y(1);
pMPF=y(2);
Ste9=y(3);
Slp1T=y(4);
Slp1=y(5);
IEP=y(6);
Rum1T=y(7);
SK=y(8);
M=y(9);
dum=y(10);

BB     = Cdc13T+Rum1T+Kdiss;
Trimer = 2.*Cdc13T*Rum1T/(BB+sqrt(BB^2-4.*Cdc13T*Rum1T));
%MPF   = (Cdc13T-pMPF)*(Cdc13T-Trimer)/Cdc13T;
MPF = Cdc13T - Trimer - pMPF*Kdiss/(Kdiss+Rum1T-Trimer);

TF     = GK(k15p*M+k15pp*SK,k16p+k16pp*MPF,J15,J16);
kwee   = kweep + (kweepp-kweep)*GK(Vawee,Viweep+Viweepp*MPF,Jawee, Jiwee);
k25    = k25p + (k25pp-k25p)*GK(Va25p+Va25pp*MPF,Vi25,Ja25,Ji25);

dy(1) 	= k1*M - (k2p+k2pp*Ste9+k2ppp*Slp1)*Cdc13T;
dy(2)  	= kwee*(Cdc13T-pMPF) - k25*pMPF - (k2p+k2pp*Ste9+k2ppp*Slp1)*pMPF;
dy(3)  	= (k3p+k3pp*Slp1)*(1.-Ste9)/(J3+1.-Ste9) -(k4p*SK+k4*MPF)*Ste9/(J4+Ste9);
dy(4)	= k5p + k5pp*(MPF^4)/(J5^4+MPF^4) - k6*Slp1T;
dy(5)	= k7*IEP*(Slp1T-Slp1)/(J7+Slp1T-Slp1) - k8*Slp1/(J8+Slp1)-k6*Slp1;
dy(6)  	= k9*(1.-IEP)*MPF/(J9+1.-IEP) - k10*IEP/(J10+IEP);
dy(7)  	= k11 - (k12+k12p*SK+k12pp*MPF)*Rum1T;
dy(8)  	= k13p + k13pp*TF - k14*SK;
dy(9)   = mu*M;
dy(10)	= K*(MPF-dum);

end

%%

function f=GK(a,b,c,d)
dum = b-a+b*c+a*d;
f=2*a*d/(dum+sqrt(dum^2-4.*(b-a)*a*d));
end

