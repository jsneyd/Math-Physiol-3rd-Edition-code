function main

clear all

global lam Vm tt n

ys=0.8;
Vm=1;
tt=[10];
n=3;
lam=ys*fb(ys);

sol=dde23(@rhs,tt,@yinit,[0 200]);

plot(sol.x,sol.y)
xlabel('time')
ylabel('CO2 concentration')

 

function dy=rhs(t,y,Z)
global lam Vm tt n
dy=[0];
dy(1) = lam - y(1)*fb(Z(1));



function out=fb(dum)
global lam Vm tt n
out = Vm*dum^n/(1+dum^n);

function out=yinit(t)
out=0.5;