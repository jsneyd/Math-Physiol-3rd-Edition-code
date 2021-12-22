function dy=ccgen5odes(t,y)

global k1 k2p k2pp k2ppp k3p k3pp k4p k4 k5p k5pp k6 k7 k8
global J3 J4 J5 J7 J8 k9 k10 mu

dy=zeros(6,1);

cycB=y(1);
cdh=y(2);
cdcT=y(3);
cdcA=y(4);
IEP=y(5);
m=y(6);

dy(1)=k1*m-(k2p+k2pp*cdh)*cycB;
dy(2)=(k3p+k3pp*cdcA)*(1.0-cdh)/(J3+1-cdh) - k4*cycB*cdh/(J4+cdh);
dy(3)=k5p+k5pp*((cycB)^4)/(J5^4+(cycB)^4)-k6*cdcT;
dy(4)=k7*IEP*(cdcT-cdcA)/(J7+cdcT-cdcA) - k8*cdcA/(J8+cdcA) - k6*cdcA;
dy(5)=k9*cycB*(1.0-IEP) - k10*IEP;
dy(6)=mu*m;
