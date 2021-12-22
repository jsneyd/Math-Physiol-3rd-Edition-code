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

% CycBT=y(1);
% pMPF=y(2);
% Cdh1=y(3);
% Cdc20T=y(4);
% Cdc20A=y(5);
% IE=y(6);
% CKIT=y(7);
% SK=y(8);
% M=y(9);
% dum=y(10);

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
MPF   = (Cdc13T-pMPF)*(Cdc13T-Trimer)/Cdc13T;

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