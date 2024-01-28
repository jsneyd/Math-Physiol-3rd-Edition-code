% this is an attempt to solve the nephron equations for loop of Henle
clear all; close all; clc;
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
% parameters
p.P = 0.9;
p.DPd = 0.15;
p.DPc = 0.22;
p.Hd = 0.1;
p.rd = 0.15;
 

p.N=51; % number of grid points
y = linspace(0,1,p.N);
 

onebyrclist = [0,0.5];
 
% make up some initial guess
Qa = -.5;
Qd = -Qa*y +(1-y);
Sd=-(Qd -1)*(p.rd*p.Hd)+1-p.DPd*p.Hd*y ; 
 
 X0=[Sd, Qa ];
 
[U]=fsolve(@(x)des(x,p),X0); 
Sd = U(1:p.N);
 
Qa = U(p.N+1);
 
tm=(1-Sd-p.DPd*p.Hd*y)/(p.rd*p.Hd);
Qd = 1+tm;

figure(1)
 
plot(y,Qd,'--')
hold on 
 
xlabel('y')
ylabel('Relative Flux')
 
box off
  
Cd=Sd./Qd;
 Ca=Cd(p.N)-p.P*(y-1)/Qa;
 
figure(2)
 
plot(y,Ca,y,Cd)
text(0.09,1.65,'C_d','fontsize',18)
text(0.18,1 ,'C_a','fontsize',18)
 
xlabel('y')
ylabel('Relative Concentration')
box off

 

%now loop on P for no ADH case:
% make up some initial guess
Qa = -.5;
Qd = -Qa*y +(1-y);
Sd=-(Qd -1)*(p.rd*p.Hd)+1-p.DPd*p.Hd*y ; 
 
 
 X0=[Sd, Qa];
Plist = [0.01:.01:1];
  
for j = 1:length(Plist)
    p.P = Plist(j);
 
[U]=fsolve(@(x)des(x,p),X0); 
% use the solution as the initial guess for the next try
X0=U;
Sd = U(1:p.N);
tm = (1-Sd-p.DPd*p.Hd*y)/(p.rd*p.Hd);
Qd = 1+ tm;
 Qa=U(p.N+1);
Cd=Sd./Qd;
 
 Ca=Cd(p.N)-p.P*(y-1)/Qa;
QdP(j) = Qd(p.N);
 Ca0(j) = Ca(1);

end
figure(4)
plot(Plist,Ca0,Plist,QdP)
xlabel('P')
axis([0 1 0 1.2])
text(0.8,.8,'C_a(0)','fontsize',18)
text(0.6,0.5,'Q_a(1)','fontsize',18)

function out = des(U,p)
 
y = linspace(0,1,p.N);
dy = y(2)-y(1);

Sd = U(1:p.N);
 
Qa = U(p.N+1);
 

Ss = p.P*(y-1) + Sd(p.N)-Sd;
tm = (1-Sd-p.DPd*p.Hd*y)/(p.rd*p.Hd);
Qs = -1-Qa  - tm ;
Qd = 1+ tm;
Fd = Ss./Qs-Sd./Qd;
 

Cd=(1+p.rd*p.Hd*(1-Qd)-p.DPd*p.Hd*y)./Qd;
Cs = (p.P+p.DPd*p.Hd).*(1-y)./(Qd+Qa)-p.rd*p.Hd;


eqSd =[Sd(1)-1,(Sd(2:p.N-1)-Sd(1:p.N-2))/(dy)-p.Hd*Fd(1:p.N-2), ...
    (Sd(p.N)-Sd( p.N-1))/(dy)-p.Hd*Fd(p.N-1)];  
% 
 
eqQa = Sd(p.N)-1-(Qa+1)*p.rd*p.Hd+p.DPd*p.Hd;
 
 
out = [eqSd, eqQa ];

end

