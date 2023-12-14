clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0);


% Hemoglobin saturation curve data
A =[0	0;
2.20994	3.07692;
3.59116	4.61538;
6.07735	6.76923;
10.49724	10.15385;
14.0884	12.3077;
19.337	15.3846;
28.453	18.7692;
40.3315	22.7692;
50	25.8462;
60.4972	30.1538;
69.8895	36;
80.1105	45.2308;
83.97790000000001	51.6923;
88.9503	61.8462;
93.3702	75.38460000000001;
95.85639999999999	87.07689999999999;
98.0663	110.462];

%MWC saturation curve
O2 = [0:.1:120]'; % O2 partial pressure
K = 26;
Y = O2.^4./(K^4+O2.^4);
K1 = 45.9;
K2 = 23.9;
K3 =243.1;
K4 = 1.52;
alp(1) = 1;
alp(2) = alp(1)/K1;
alp(3)=alp(2)/K2;
alp(4) = alp(3)/K3;
alp(5) = alp(4)/K4;
vec = [ones(length(O2),1),O2,O2.^2,O2.^3,O2.^4];

ints = [0:4];
intalps = ints.*alp;
SatY = (vec*intalps')./(vec*alp');

 
figure(1)
plot(A(:,2),A(:,1),'*',O2,100*Y,O2,25*SatY,'--')
xlabel('Oxygen Partial Pressure (mm Hg)')
ylabel('Saturation (%)')
box off

figure(2)

x = [0:.01:10];
Y=(x+3*x.^2+3*x.^3+x.^4)./(1+4*x+6*x.^2+4*x.^3+x.^4);

K1 = 1000;
K2 = 1250;
K3 = 4.54;
YMWC = (O2/K1.*(1+O2/K1).^3 + O2./K3.*(1+O2./K3).^3./K2)./((1+O2/K1).^4 +(1+O2./K3).^4./K2);

plot(A(:,2),A(:,1),'*',O2,25*SatY,'--',O2,100*YMWC,'b')
xlabel('Oxygen Partial Pressure (mm Hg)')
ylabel('Saturation (%)')
box off
