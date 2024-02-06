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
p.Ka = 0.2;

p.N=51; % number of grid points
p.y = linspace(0,1,p.N);
 

onebyrclist = [0,0.5];
 
% make up some initial guess
Qa = -.5;
Qd = -Qa*p.y +(1-p.y);
Sd=-(Qd -1)*(p.rd*p.Hd)+1-p.DPd*p.Hd*p.y ; 
 
 X0=[Sd, Qa ];
 
[U]=fsolve(@(x)des(x,p),X0); 
Sd = U(1:p.N);
 
Qa = U(p.N+1);
 
Sd = U(1:p.N);
Qa = U(p.N+1);
[Qd,Qs,Sd,Sa,Ss,Cd,Ca,Cs]=get_fs(U, p);
 
 
figure(1)
 
plot(p.y,Qd,'--')
hold on 
 
xlabel('y')
ylabel('Relative Flux, Q_d')
 
box off
 
figure(2)
 
plot(p.y,Ca,p.y,Cd)
 
legend('boxoff')
legend('C_a','C_d','location','northwest') 
xlabel('y')
ylabel('Relative Concentration')
box off
 
tic
%now loop on P for no ADH case:
% make up an initial guess
Qa = -.5;
Qd = -Qa*p.y +(1-p.y);
Sd=-(Qd -1)*(p.rd*p.Hd)+1-p.DPd*p.Hd*p.y ; 
 
 
 X0=[Sd, Qa];
Plist = [0.01:.01:1];
  
for j = 1:length(Plist)
    p.P = Plist(j);
 
[U]=fsolve(@(x)des(x,p),X0); 
% use the solution as the initial guess for the next try
X0=U;
Sd = U(1:p.N);
Qa = U(p.N+1);
[Qd,Qs,Sd,Sa,Ss,Cd,Ca,Cs]=get_fs(U, p);
 
 
QdP(j) = Qd(p.N);
 Ca0(j) = Ca(1);
Cd1(j) = Cd(p.N);
end
toc
figure(3)
plot(Plist,Ca0,Plist,QdP, Plist,Cd1)
xlabel('P')
%axis([0 1 0 1.3])
%text(0.8,.8,'C_a(0)','fontsize',18)
%text(0.6,0.5,'Q_d(1)','fontsize',18)
legend('boxoff')

legend('C_a(0)','Q_d(1)','C_d(1)','location','northwest')

function out = des(U,p)
 
 
dy = p.y(2)-p.y(1);
Sd = U(1:p.N);
Qa = U(p.N+1);
[Qd,Qs,Sd,Sa,Ss,Cd,Ca,Cs]=get_fs(U ,p);
 
 
Fd = Ss./Qs-Sd./Qd;
 
eqSd =[Sd(1)-1,(Sd(2:p.N-1)-Sd(1:p.N-2))/(dy)-p.Hd*Fd(1:p.N-2), ...
    (Sd(p.N)-Sd( p.N-1))/(dy)-p.Hd*Fd(p.N-1)];  
%  
eqQa = Sd(p.N)-1-(Qa+1)*p.rd*p.Hd+p.DPd*p.Hd;

out = [eqSd, eqQa ];

end

function [Qd,Qs,Sd,Sa,Ss,Cd,Ca,Cs]= get_fs(U, p)
Sd = U(1:p.N);
Qa = U(p.N+1);
 
tm = (1-Sd-p.DPd*p.Hd*p.y)/(p.rd*p.Hd);
Qd = 1+ tm;
Qs=-Qd-Qa;
Cd=Sd./Qd;
Ca=Cd(p.N)-p.P*(p.y-1)/Qa;
 
Sa =Qa*Ca;
Ss = -Sa -Sd;
Cs=Ss./Qs;
 
end
 