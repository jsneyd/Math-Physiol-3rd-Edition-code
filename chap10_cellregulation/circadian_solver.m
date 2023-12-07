%this routine solves a system of odes;  Here is the Tyson et al. circadian
%model
function circadian_clock
clear
global vm vp km k3 k1 k2 J K A

%specify parameter values
vm = 1;
vp = 0.5;
km = .1;
k3 = .1;
k1 = 10;
k2 = 0.03;
J = 0.05;
K = 200;
A = 0.1;  %Pcrit

P = [ 0.001:.001:20];
q = 2./(1+sqrt(1+8*K*P));
P1 = q.*P;
P2 = ((1-q)/2).*P;

M1  = vm./((1+(P2/A).^2) * km);

M2 =((k1*P1+2*k2*P2)./(J+P)+k3*P)/vp;

%specify the output points
tspan = [-24:.1:50];

%specify initial data

y0 = [0,0];



[T,Y] = ode23(@deRHS,tspan, y0);  
figure(1)
plot(T,Y(:,1),'r',T,Y(:,2),'b','linewidth',2)
xlabel('time (hours)','fontsize',16)
text(16,3,'PER','fontsize',16)
text(11,1.5,'mRNA','fontsize',16)
axis([0 50 0 3.5])
%legend('mRNA','Per')


figure(2)
loglog(P,M1,'b--',P,M2,'r--',Y(:,2),Y(:,1),'linewidth',2)
ylabel('mRNA','fontsize',16)
xlabel('Per','fontsize',16)
text(2.5,.1,'dM/dt = 0','fontsize',16)
text(2,4,'dP/dt = 0','fontsize',16)
axis([.001 100 .01 100])

function F_prime=deRHS(t,y)
global vm vp km k3 k1 k2 J K A

M = y(1);
P = y(2);
q = 2/(1+sqrt(1+8*K*P));
P1 = q*P;
P2 = ((1-q)/2)*P;

Xp = vm/(1+(P2/A)^2) - km*M;

Yp = vp*M-(k1*P1+2*k2*P2)/(J+P)-k3*P;

F_prime = [Xp;Yp];
