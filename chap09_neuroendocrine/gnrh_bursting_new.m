% Bursting in GnRH neurons.
% Simplified version of the model of Duan et al, Journal of Theoretical Biology 276 (2011) 22–34
% For Keener and Sneyd, third edition

%this version uses the Li-Rinzel model for calcium relese rather than a
%modal model.


close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 

%parameters
clear all; close all; clc
format long;

global p

init = [     
 -47.439744033440896
   0.150875059704950
   0.068499547036174
   0.042896908631945
   10
   0.061113223260552
   0.000000002523211
   0.000049694789395
];

p.IP3=0.15;  %IP3 level is a parameter
p.Kf=1.92*10^(-1);
p.Kf=1.92*10^(-2);
 

% Ca and IPR parameters
p.tau_max=1000; p.Ktau=0.1;
p.tauP=0.5;  p.Kc=0.2; p.Kh=0.08;
p.tauh=500;

p.gnaf= 280; p.gsk=0.3; 
p.gw=950;  p.gcal=0.05; 
p.Prate=1; p.Vp=0.0042; 

p.gkm=8; p.gleak=0.04; 
p.Cm=16; p.Vna=60; p.Vca=100; p.Vk=-80; p.Vleak=100; %Resting potentials  
p.k_sk=1; p.n_sk=2; p.k11=1*10^(-7); p.k_11=1.2; p.k22=0.5; p.k33=3*10^(-5); %Isk and Iw channel parameters

p.rho=0.5; p.gamma=27;%Compartment converting parameters and sometimes rho values can be modified for TTX result
p.a1=1*10^(-4); p.a2=35; p.a3=300; p.a4=7; p.a5=35; %SERCA pump parameters 
p.alpha=4.8*10^(-3);  p.Kp=0.425; p.beta=2*10^(-5); p.Jer=4*10^(-8); p.Knaca=0.05; %Fluxes parameters

% from Keizer-DeYoung
p.k1 = 400;
p.k2 = 0.2;
p.k3 = 400;
p.k4 = 0.2;
p.k5 = 20;
p.km1 = 52;
p.km2 = 0.21;
p.km3 = 377.2;
p.km4 = 0.0289;
p.km5 = 1.64;
p.K1=p.km1/p.k1;
p.K5=p.km5/p.k5;
p.K2 = p.km2/p.k2;
p.K3=p.km3/p.k3;
p.K4=p.km4/p.k4;
p.fer=0.01;
  p.fcyt=0.01;

tic
[Time,soln1]=ode15s(@Model,[0 60000],init); 
toc
V = soln1(:,1);
Hnaf = soln1(:,2);
Nkm = soln1(:,3);
c = soln1(:,4);
ce = soln1(:,5);
h = soln1(:,6);
Ow = soln1(:,7);
Ow_star = soln1(:,8);

X = Ow + Ow_star;

tsec = Time/1000;
% figure(10)
% subplot(3,1,1);
% plot(tsec,V); ylabel('V')  ;
% subplot(3,1,2); plot(tsec, c); ylabel('c');
% subplot(3,1,3); plot(tsec,ce); xlabel('time (s)'); ylabel('ce');
% 
figure(1)
%
plot(tsec,V)
ylabel('V (mV')
xlabel('t (s)')
axis([39 46 -65 30])

yyaxis('right')
plot(tsec,c)
ylabel('[Ca^{++}] (\muM')

figure(11)
plot(c,V)

% Create a table with the data and variable names
T = table(tsec,V,c,ce,X, 'VariableNames', {'tsec', 'V', 'c', 'ce','X'});
% Write data to text file
writetable(T, 'temp.txt')

%%
function sys = Model(Time,S,p)
global p

V = S(1); Hnaf = S(2); Nkm = S(3); 
c = S(4); ce = S(5); h = S(6); Ow = S(7); Ow_star = S(8);

% Current submodel
Mnaf_inf = 1/(1+exp(-(V+40.0)/4.3));
Hnaf_inf = 1/(1+exp((V+66.1)/10.8));
Nkm_inf = 1/(1+exp(-(V+37)/4));
Mcal_inf = 1/(1+exp(-(V+30)/2));

Thnaf = 75/(exp((V+80)/19)+2*exp(-2*(V+80)/19));
Tnkm = 11.5/(exp((V+30)/15)+exp(-(V+30)/15));

Inaf = p.gnaf*Mnaf_inf^3*Hnaf*(V-p.Vna);
Ical = p.gcal*Mcal_inf^2*(V-p.Vca);
Ileak = p.gleak*(V-p.Vleak);

Ikm = p.gkm*Nkm*(V-p.Vk);
Isk = p.gsk*(c^(p.n_sk)/(c^(p.n_sk)+p.k_sk^(p.n_sk)))*(V-p.Vk);
Iw = p.gw*(Ow+Ow_star)*(V-p.Vk);

Iionic = Inaf+Ikm+Ical+Ileak;
Total_calcium = Ical;

% calcium submodel
Jserca = p.Prate*(c-p.a1*ce)/(p.a2+p.a3*c+p.a4*ce+p.a5*c*ce);
Jin = -p.alpha*(Total_calcium)+p.beta*p.IP3;
Jpm = p.Vp*c^2/(p.Kp^2+c^2);

% this is a basic Li-Rinzel model (Keizer-deyoung simplified)
Po = (p.IP3*c*h./((p.IP3+p.K1).*(c+p.K5))).^3; %open probability
   
Jrel = (p.Kf*Po+p.Jer)*(ce-c); 
 
dVdt = (-1/p.Cm)*(Iionic+Isk+Iw);
dHnafdt = (Hnaf_inf-Hnaf)/Thnaf;
dNkmdt = (Nkm_inf-Nkm)/Tnkm;
 

dcdt = p.rho*(Jin - Jpm) + (Jrel  - Jserca); 
%dhdt  = 0.001*((h_inf-h)/tau);  %old version
dhdt=(p.K2-(c+p.K2)*h)/p.tauh;
 %new verson
dcedt = p.gamma*(Jserca - Jrel);
dOwdt = -Ow*p.k_11-Ow*p.k22+p.k11*c*(1-Ow-Ow_star); 
 
dOw_stardt = p.k22*Ow-p.k33*Ow_star;

sys = [dVdt;...
    dHnafdt;dNkmdt;...
    dcdt;dcedt;dhdt;...
    dOwdt;dOw_stardt...
    ];
end
