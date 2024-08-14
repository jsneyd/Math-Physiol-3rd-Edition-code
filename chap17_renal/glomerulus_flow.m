% glomerular filtration, unregulated
clear all
close all
clc

 
global  KfL pii Ra Re Qi Pe P2 Pa
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
  
%parameters (typical values)
P1 = 60;
P2 = 18;
Pa = 100;
Pe = 18;
Pd = 18;% = P2 so that Rd = 0
pii = 25; %mm Hg
Qi = 650;
Qd = 125;
Qe0 = Qi-Qd
alpha = pii/(P1-P2)
%use typical parameters to find Ra, Re, Rd and KfL
Ra = (Pa-P1)/Qi
Re = (P1-Pe)/Qe0
Rd = (P2-Pd)/Qd
Qe=Qe0;
KfL = -(Qe./Qi +alpha*log((Qe./Qi-alpha)./(1-alpha)) -1) *(alpha *Qi)/pii

Palist = [20:1:160];  
for j = 1: length(Palist)
   Pa = Palist(j);

   qemx = (Pa-Pe)/(Ra+Re);
   qemin = (P2*Ra - Pe*Ra - Re*pii + sqrt((P2^2  - 2*P2*Pe  + Pe^2)*Ra^2  ...
       - 2*(P2 - 2*Pa +Pe)*Re*Ra*pii + Re^2*pii^2)  )/(2*Re*Ra)+0.1;
  
   Qe(j)=bisect(@F,qemin,qemx); 
 
end

 figure(1)

 Qi= (Palist-Pe-Re*Qe)/Ra;
 Qd=Qi-Qe;
 ndx=find(Qd>0);
    plot( Palist(ndx),Qd(ndx))
    xlabel('Arterial pressure, P_a (mm Hg)')
    ylabel('Glomerular filtration rate, Q_d (ml/min)')
    yyaxis right 
    plot(Palist(ndx),Qe(ndx))
    ylabel('Efferent arteriole flow rate, Q_e (ml/min)')

   % axis([0 150 0 1000])
% figure(2)
%     plot(Pa,Qdlist,  Pa,0.2*Qi)
%     xlabel('arterial pressure')
%     legend('Q_d','Q_i')
%     axis([20 120 0 150])

%%
function out = F(Qe)
    global  KfL pii Ra Re P2 Pe Pa
    Qi= (Pa-Pe-Re*Qe)/Ra;
    P1 = -Ra*Qi+Pa;
    alpha = pii/(P1-P2);
     
    out = Qe +alpha.*Qi.*log((Qe./Qi-alpha)./(1-alpha)) -Qi + KfL*pii./(alpha);
   % out=  ((Qe./Qi-alpha) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% the bisection algorithm
function root = bisect(feval,a,b)

ul = a;
fl = feval(ul);
uu = b;
fu = feval(uu);

% we make the assumption, without checking, that
% fu*fl<0

% if not, the algoritm fails to find a root.

N = 25;  % number of iterates
% the main bisection algorithm
for j = 1:N
u = (ul+uu)/2;
fc = feval(u);
ftest = (fc*fl>=0);
ul = ftest*u+(1-ftest)*ul;
fl = ftest*fc + (1-ftest)*fl;

uu = (1-ftest)*u+ ftest*uu;
fu = (1-ftest)*fc + ftest*fu;
end
root = u;
end


 