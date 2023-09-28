%model of l-type calcium channel

function L_type_Ca

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0); 


%parameters
Vm = -75;
c = 0;
p.alpha=2*exp(0.012*(Vm-35)); 
p.beta = 0.0882*exp(-0.05*(Vm-35));
p.gamma = 0.44*c;
p.omega = 0.02;
p.f=0.85;
p.g=2;

 Vmlist=[-65,-55,-45,-35,-25,-15,-5];
 for j=1:length(Vmlist)
     Vm=Vmlist(j);
     p.alpha=2*exp(0.012*(Vm-35)); 
p.beta = 0.0882*exp(-0.05*(Vm-35));
dt = 0.01;
tend=15;

init = [0.965,0.0248,0,0,0.0102,0];
tspan = [0:dt:tend];
[t,sol] = ode15s(@(t,x)Ltrhs(t,x,p),tspan,init);


  op=sol(:,5);
  c = sol(:,6)
  figure(1)

  plot(t,op)
  hold on


  figure(2)
  plot(t,c)
  hold on

 end

 


    
    
    
function out=Ltrhs(t,x,p)
 
n3=x(1);
n4=x(2);
c3=x(3);
c4=x(4);
op=x(5);
c = x(6);
p.gamma = 0.44*c;

n3p=-(p.alpha+p.gamma)*n3+p.beta*n4+p.omega*c3;
n4p=-(p.f+p.beta+2*p.gamma)*n4+p.alpha*n3+p.g*op+p.omega/2*c4;
c3p=-(2*p.alpha+p.omega)*c3+p.gamma*n3+p.beta/2*c4; 
c4p=-(p.beta/2+p.omega)*c4+2*p.gamma*n4+2*p.alpha*c3;
opp = -p.g*op+p.f*n4;
cp = op-c;


out=[n3p;n4p;c3p;c4p;opp;cp];

 