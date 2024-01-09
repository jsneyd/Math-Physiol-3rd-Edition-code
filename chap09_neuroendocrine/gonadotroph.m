%-----------------------------------
%
% Pituitary cell model using a common core of ionic currents.
% Adding ER calcium, a cytosolic calcium subspace, and IP3 receptor dynamics.
% Adapted from the original code (Fletcher, Bertram, Stojilkovic, Mol Cell Endocrin 463 (2018))
% for Keener and Sneyd, Math Physiology, third edition.
%
% ----------------------------------

close all
clear all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);

global ip3b ip3p tpulse gnav gca gk gkca gL gkir ek eca ena vL alpha Cm sigmaER  sigmaD kp tau 
global tauhna taun kca fcyt Kdiff  kleak kf ki  ka Kd  Aip3 fer Vpmca  Kpmca Vserca Kserca pnorm

% parameters
% background IP3 and pulse IP3. Vary  these to get different bursting patterns.
ip3b=0;
ip3mx=[0.1,0.4,0.7,2];

tpulse = 5000;


gnav=20;
gca=1.5;
gk=5;
gkca=2.5;
gL=0.1;
gkir=0.9;

ek=-75;
eca=60;
ena=75;
vL=0;
alpha=0.0015;
Cm=6;
sigmaER=39;
sigmaD=39;
kp=0.04;
tauhna=2;
taun=20;
tau=175000;
pnorm=4*tau^2*exp(-2);  % sets a time scale for the ip3 input 
Vserca=0.16;
Kserca=0.2;
Kdiff=0.6;
kca=0.4;
fer=0.01;
fcyt=0.01;
Vpmca=0.02;
Kpmca=0.1;
kleak=0.0002;
kf=20;
ki=1;
ka=.8;
Kd=.8;
Aip3=0.002;
 
for j = 1:4
    ip3p=ip3mx(j);; % This sets the maximum value of the IP3 pulse input
    %%%%%  initial Conditions %%%%%%%%%%%%%%%
    v0=-60.58;
    c0=0.0881;
    cer0=124.49;
    cd0=0.130;
    
    hna0=0.637;
    n0=0.0;
    
    hip3r0=0.862;
    
    init =[v0,c0,cer0,cd0,hna0,n0,hip3r0];
    
    total=60000;
    tstep =  1;
    tic
    %specify the output points
    tspan = [0:tstep:total];
    [T,S] = ode23(@deRHS,tspan,init  , odeset('maxstep',10));  
    toc
    figure(2*j-1)
    
    plot(T/1000,S(:,1))
    ylabel('V (mV)')
    xlabel('t (s)')
    yyaxis('right')
    plot(T/1000, ip3b+ip3p*ip3pulse(T))
    ylabel('IP_3')
    
    figure(2*j)
    plot(T/1000,S(:,2))
    xlabel('t (s)')
    ylabel('c(\muM}')
    yyaxis('right')
    plot(T/1000,S(:,3))
    ylabel('C_e (\muM)')

end


%% The equations
function s_prime=deRHS(t,sol) 

global ip3b ip3p  gnav gca gk gkca gL gkir ek eca ena vL alpha   Cm  sigmaER  sigmaD 
global tauhna taun kca fcyt Kdiff  kleak kf ki  ka Kd  Aip3 fer Vpmca  Kpmca Vserca Kserca

%there are 7 variables
v=sol(1);
c=sol(2);
cer = sol(3);
cd = sol(4);
hna = sol(5);
n = sol(6);
hip3r = sol(7);

vKdrive=v-ek;
vnadrive=v-ena;
vcadrive=v-eca;
vnsdrive=v-vL;

%%%%%% cell geometry %%%%%%%%%%%%

%%%%%%% smooth pulse profiles %%%%%%%%%%%%
% time to peak is 2*tau, peak value is 1. 


%%%%%% sodium %%%%%%%%%%%%%%%%%%%%%%%%%

% ina -- ttx-sensitive

mnainf = 1/(1+exp((-15-v)/5));
hnainf = 1/(1+exp((55+v)/10));
in = gnav*mnainf^3*hna*vnadrive;

%%%%%% calcium %%%%%%%%%%%%%%%%%%%%%%%%
% ica
mcainf = 1/(1+exp((-15-v)/12));
ica = gca*mcainf*vcadrive;

%%%%%% potassium %%%%%%%%%%%%%%%%%%%%%%
% ik

ninf = 1/(1+exp((-v)/5));
iKdr = gk*n*vKdrive;

% ikca -- sk, ik, and/or bk far from ca2+ channels

c2=c^2;

nkcainf=c2/(c2+kca^2);
ikca = gkca*nkcainf*vKdrive;

% ikir and GIRK
nkirinf = 1/(1+exp((55+v)/10));
ikir = gkir*nkirinf*vKdrive;

ik = ikca + iKdr  + ikir;

%%%%%%%% non-specific leak %%%%%%%%%%%%%%%%%

iL=gL*(v-vL);

%%%%%% calcium handling %%%%%%%%%%%%%%%

jin = -alpha*ica;

% PMCA
jpmca=Vpmca*c2/(c2+Kpmca^2) ;

% SERCA
jserca=Vserca*c2/(c2+Kserca^2);

%jrel consists of ER leak and IP3 receptor

% background IP3 + pulse
ip3=ip3b+ip3p*ip3pulse(t);


% IPR stuff
num=(ip3*cd*hip3r);
num3= num^3;
denom=(ip3+ki)*(cd+ka);
denom3=denom^3;
jip3r=kf*num3/denom3*(cer - cd);

jleak=kleak*(cer - cd);
jrel = jleak+jip3r;

% calcium subspace
jdiff=Kdiff*(cd-c);

%%%%% Equations %%%%%%%%%%%%%%%%%%%%%%%
vp= -(in+ica+ik+iL)/Cm;
cp = fcyt*(jin-jpmca+jdiff-jserca);
cerp = fer* sigmaER*(jserca-jrel);
cdp = fcyt* sigmaD*(jrel-jdiff);
hnap= (hnainf-hna)/tauhna;
np= (ninf-n)/taun;
hip3rp = Aip3*(Kd-(cd + Kd)*hip3r);

s_prime=[vp;cp;cerp;cdp;hnap;np;hip3rp];
end


%% IP3 pulse
function ip3out=ip3pulse(t)
global tpulse   pnorm kp   

delpulse=t-tpulse;
pulse=delpulse.^2./pnorm.*(delpulse>0);
ip3pulse=pulse./(pulse+kp)*(1+kp);
ip3out = ip3pulse;
end