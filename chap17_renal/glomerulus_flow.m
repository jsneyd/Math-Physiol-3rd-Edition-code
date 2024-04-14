% glomerular filtration, unregulated
global alpha KfL Qd pii
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
  
%parameters (typical values)
P1 = 60;
P2 = 18;
Pa = 100;
Pe = 0;
Pd = 18;
pii = 25; %mm Hg
Qi = 650;
Qd = 125;
Qe = Qi-Qd;
alpha = pii/(P1-P2);
%use typical parameters to find Ra, Re, Rd and KfL
Ra = (Pa-P1)/Qi;
Re = (P1-Pe)/Qe;
Rd = (P2-Pd)/Qd;

KfL = -(Qe./Qi +alpha*log((Qe./Qi-alpha)./(1-alpha)) -1) *(alpha *Qi)/pii;
% require Qi>Qd/(1-alpha)

Qdlist = [25:1:150];  
for j = 1:length(Qdlist)
    Qd = Qdlist(j);
    Qi(j)=fsolve(@F,100);
end
Qe = Qi-Qdlist;
P1 = Re*Qe+Pe
Pa=P1+Ra*Qi;
P2 = Pd+Rd*Qdlist

figure(1)
    plot(Qdlist,Qi,Qdlist,Qe)
    xlabel('Q_d')
    legend('Q_i','Pa')
    axis([0 150 0 1000])
figure(2)
    plot(Pa,Qdlist,  Pa,0.2*Qi)
    xlabel('arterial pressure')
    legend('Q_d','Q_i')
    axis([20 120 0 150])

%%
function out = F(Qi)
    global alpha  Qd KfL pii
    Qe = Qi-Qd;
    out = Qe +alpha.*Qi.*log((Qe./Qi-alpha)./(1-alpha)) -Qi + KfL*pii./(alpha);
end

 