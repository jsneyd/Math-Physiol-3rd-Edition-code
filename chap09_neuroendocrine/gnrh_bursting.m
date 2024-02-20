
%  -------------------------------------------------------------------
%
%   Solve the model of bursting in GnRH neurons.
%   Simplified version of the model of Duan et al, Journal of Theoretical Biology 276 (2011) 22–34
%
%   For Chapter 9, Section 9.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 

% parameters
global p
p.IP3=0.15;
p.Kf=1;

% Ca and IPR parameters
p.tau_max=10; p.Ktau=0.1;
p.tauP=0.5;  p.Kc=0.2; p.Kh=0.08; p.Kp=0.425;
p.delta=0.5; p.gamma=27;
p.alpha0=4.8*10^(-3);   p.alpha1=2*10^(-5); 
p.Vp=0.0042;
p.Vserca=0.1;
p.kserca=0.2;
p.Kbar=0.00001957;

% conductances
p.gnaf= 280; p.gsk=0.3; 
p.gw=950;  p.gcal=0.05;   
p.gkm=8; p.gleak=0.04; 

% Resting potentials
p.Cm=16; p.Vna=60; p.Vca=100; 
p.Vk=-80; p.Vleak=100;   

% Isk and Iw channel parameters
p.k_sk=1; p.n_sk=2; p.k11=1*10^(-7); 
p.k_11=1.2; p.k22=0.5; p.k33=3*10^(-5); 



% Solve the ODEs
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
[Time,soln1]=ode15s(@Model,[0 60000],init); 

V = soln1(:,1);
Hnaf = soln1(:,2);
Nkm = soln1(:,3);
c = soln1(:,4);
ce = soln1(:,5);
h = soln1(:,6);
Ow = soln1(:,7);
Ow_star = soln1(:,8);

X = Ow + Ow_star;

figure(1)
%
plot(Time/1000,V)
ylabel('V (mV')
xlabel('t (s)')
axis([52 58 -65 30])

yyaxis('right')
plot(Time/1000,c)
ylabel('[Ca^{++}] (\muM')

% figure(2)
% plot(c,ce)
% ylabel('c_e (mM')
% xlabel('c (mM')
% 
% figure(3)
% plot(c,V)
% ylabel('V (mV')
% xlabel('c (mM')
% 
% figure(4)
% plot(tsec,c,tsec,ce)
% legend('boxoff')
% legend('c','c_e')
% xlabel('t (s)')

%% stuff for external plotting
% Create a table with the data and variable names
% T = table(Time/1000,V,c,ce,X, 'VariableNames', {'tsec', 'V', 'c', 'ce','X'});
% Write data to text file
% writetable(T, 'temp.txt')

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

% calcium submodel
phi_c=c^4/(c^4+p.Kc^4);
phi_p=p.IP3^2/(p.Kp^2+p.IP3^2);
phi_p_down=p.Kp^2/(p.Kp^2+p.IP3^2);

h_inf=p.Kh^4/(p.Kh^4+c^4);
tau = p.tau_max*p.Ktau^4/(p.Ktau^4+c^4);
beta = phi_p*phi_c*h;
alpha = phi_p_down*(1-phi_c*h_inf);
Po=beta/(beta + 0.4*(beta+alpha));

JIPR= p.Kf*Po*(ce-c); 
Jserca = p.Vserca*(c*c-p.Kbar*ce*ce)/(c*c+p.kserca*p.kserca);
Jin = -p.alpha0*(Ical)+p.alpha1*p.IP3;
Jpm = p.Vp*c^2/(p.Kp^2+c^2);

dVdt = (-1/p.Cm)*(Iionic+Isk+Iw);
dHnafdt = (Hnaf_inf-Hnaf)/Thnaf;
dNkmdt = (Nkm_inf-Nkm)/Tnkm;
dcdt = p.delta*(Jin - Jpm) + (JIPR - Jserca); 
dhdt  = 0.001*(h_inf-h)/tau;   % Note the factor of 1/1000, since time is in ms, not s.
dcedt = p.gamma*(Jserca - JIPR);
dOwdt = -Ow*p.k_11-Ow*p.k22+p.k11*c*(1-Ow-Ow_star);
dOw_stardt = p.k22*Ow-p.k33*Ow_star;

sys = [dVdt;...
    dHnafdt;dNkmdt;...
    dcdt;dcedt;dhdt;...
    dOwdt;dOw_stardt...
    ];
end
