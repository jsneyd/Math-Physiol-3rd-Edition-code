function main
% the mackey-Glass equations, solved using dde23

clear all
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0);

global lam Vm tt n

ys=0.8;
Vm=1;
tt=[10];
n=3;
lam=ys*fb(ys);

% plot stability curve

y = [0:.01:1.3];

w = sqrt((y.*fbp(y)).^2)-(fb(y).^2);
s = (pi+atan(-w./fb(y)))./w;

figure(1)
plot(y,w, y,s)
axis([0 1.3 0 10])
text(0.65,6,'unstable','fontsize',18)
text(0.2,4,'stable','fontsize',18)
text(0.62,5,'\sigma','fontsize',20)
text(0.8,1,'\omega','fontsize',20)

sol=dde23(@rhs,tt,@yinit,[-20 100]);
 figure(2)

plot(sol.x,sol.y, sol.x,(-sol.yp+lam)./sol.y,'--')
xlabel('time')
legend('boxoff')
legend('CO_2 concentration','Ventilation rate')
axis([0 100 0 3])
 box off

function dy=rhs(t,y,Z)
global lam  
dy=[0];
dy(1) = lam - y(1)*fb(Z(1));



function out=fb(dum)
global  Vm   n
out = Vm*dum.^n./(1+dum.^n);

function out=fbp(dum)
global  Vm   n
out = Vm*n.*dum.^(n-1)./(1+dum.^n).^2;


function out=yinit(t)
out=0.5;