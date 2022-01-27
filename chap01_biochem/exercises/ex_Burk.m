set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

   
   N=7;
S=[.1;.2;.5;1;2;3.5;5];
V=[.04;.08;.16;.26;.32;.39;.42];
A=[ones(1,N)',-V./S];
B=transpose(A)*A;
c=transpose(A)*V;
b=B\c;
s=[0:100]/20;
g=b(1)*s./(s+b(2));
figure(1);
plot(S,V,'k*',s,g,'r','linewidth',2)
xlabel('Substrate Concentration S','fontsize',16)
ylabel('Reaction Velocity V','fontsize',16)
t=-.004*[0:100];
f=b(1)+t*b(2);

figure(2)
plot(t,f,'g','linewidth',2)
hold on
plot(-V./S,V,'k*')
hold off
xlabel('V/S','fontsize',16)
ylabel('V','fontsize',16)
