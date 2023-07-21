% plot of a_c,min(K_b)
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 

Kb  = [0.01:.01:10];

ac = -((3*Kb.^2+2*Kb).*log((Kb+1)./Kb)-(3*Kb+1/2))./((Kb+1/2).*log((Kb+1)./Kb)-1)/2;

figure(1)
semilogx(Kb,ac,'r')
xlabel('K_b')
ylabel('a_{c,min}')
axis([0.01 10 0.1 0.5])