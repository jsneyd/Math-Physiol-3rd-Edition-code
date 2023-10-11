% reentry loop map

 clear all
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)

L=1;

Tr = 3/4;


t = [0:.01:5];

dt = [Tr:0.01:10];
c = 1/2+2*(dt-Tr)./(1+(dt-Tr));
Lbyc = L./c;

ndx = min(find(Lbyc<=Tr));
figure(1)
T=1.75;

plot(dt(1:ndx),Lbyc(1:ndx),dt(ndx:end),Lbyc(ndx:end),':',t,t,'--',t,Tr*ones(length(t),1),':','linewidth',3)
hold on
plot([0,dt(ndx)],[T,T], ':',[dt(ndx),5],[T,T],'linewidth',3)
axis([0 2.5 0 2.5])
xlabel('\Delta T_n')
ylabel('\Delta T_{n+1}')
text(0.2,0.9,'T_r','fontsize',20)
text(0.2,1.85,'T','fontsize',20)
text(0.8,2.1,'L/c','fontsize',20)
 hold off

 T=1.25;
 figure(2)

 plot(dt(1:ndx),Lbyc(1:ndx),dt(ndx:end),Lbyc(ndx:end),':',t,t,'--',t,Tr*ones(length(t),1),':','linewidth',3)
hold on
plot([0,dt(ndx)],[T,T], ':',[dt(ndx),5],[T,T],'linewidth',3)
axis([0 2.5 0 2.5])
xlabel('\Delta T_n')
ylabel('\Delta T_{n+1}')
text(0.2,0.9,'T_r','fontsize',20)
text(0.2,1.35,'T','fontsize',20)
text(0.8,2.1,'L/c','fontsize',20)

 hold off
