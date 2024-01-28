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
 

onebyrclist0 = [0,0.5];
for j = 1:2
    p.onebyrc=onebyrclist0(j);

% make up some initial guess
Qa = -.5;
Qd = -Qa*y +(1-y);
Sd=-(Qd -1)*(p.rd*p.Hd)+1-p.DPd*p.Hd*y ; 
Qc1=0.1;
Qc = Qc1*y -Qa*(1-y);

 X0=[Sd,Qc,Qa,Qc1];
 
[U]=fsolve(@(x)des(x,p),X0); 
Sd = U(1:p.N);
Qc = U(p.N+1:2*p.N);
Qa = U(2*p.N+1);
Qc1 = U(2*p.N+2);
V = get_conc(U,p);
 
Ss = V(1:p.N);
Qs = V(p.N+1:2*p.N);
Qd = V(2*p.N+1:3*p.N);
Cd = V(3*p.N+1:4*p.N);
Cs= V(4*p.N+1:5*p.N);
Ca =  V(5*p.N+1:6*p.N);
Cc=V(6*p.N+1:7*p.N); 

figure(1)
if (j==1)
plot(y,Qd,'--')
hold on 
else 
  plot(y,Qd,y,Qc,y,-Qa*ones(p.N,1))
hold off
xlabel('y')
ylabel('Relative Flux')
text(0.6,0.75,'Q_d','fontsize',18)
text(0.6,0.58,'-Q_a','fontsize',18)
text(0.8,0.3,'Q_d','fontsize',18)
text(0.8,0.1,'Q_c','fontsize',18)
box off
end
 
if(j==1)
figure(2)
 
plot(y,Cd,y,Ca,y,Cc)
text(0.09,1.65,'C_d','fontsize',18)
text(0.18,1 ,'C_a','fontsize',18)
text(0.4,0.7,'C_c','fontsize',18)
xlabel('y')
ylabel('Relative Concentration')
box off
else
    figure(3)
   plot(y,Cd,y,Ca,y,Cc)
   text(0.4,1.45,'C_d','fontsize',18)
text(0.4,.9 ,'C_a','fontsize',18)
text(0.4,0.35,'C_c','fontsize',18)
xlabel('y')
ylabel('Relative Concentration')
box off
%axis([0 1 0 2])

end
end

%now loop on rc for no ADH case:
% make up some initial guess
Qa = -.5;
Qd = -Qa*y +(1-y);
Sd=-(Qd -1)*(p.rd*p.Hd)+1-p.DPd*p.Hd*y ; 
Qc1=0.1;
Qc = Qc1*y -Qa*(1-y);

 
 X0=[Sd,Qc,Qa,Qc1];
 
onebyrclist1 = [0.01:0.05:2];
  
for k = 1:length(onebyrclist1)
   p.onebyrc = onebyrclist1(k)
 
[U]=fsolve(@(x)des(x,p),X0); 
% use the solution as the initial guess for the next try
X0=U;
 
Sd = U(1:p.N);
Qc = U(p.N+1:2*p.N);
Qa = U(2*p.N+1);
Qc1 = U(2*p.N+2);

V = get_conc(U,p);
 
Ss = V(1:p.N);
Qs = V(p.N+1:2*p.N);
Qd = V(2*p.N+1:3*p.N);
Cd = V(3*p.N+1:4*p.N);
Cs= V(4*p.N+1:5*p.N);
Ca =  V(5*p.N+1:6*p.N);
Cc=V(6*p.N+1:7*p.N); 
 
QdP(k) = Qd(p.N);

 Qc1j(k) = Qc1;
 Cd1(k) = Sd(p.N)/Qd(p.N);
 
 Cc1(k)=Cc(p.N);
 Cc0(k)=Cc(1);


end
figure(4)
plot(onebyrclist1,Qc1j,onebyrclist1,QdP)
xlabel('1/\rho_c')
ylabel('Flow Rate')
text(0.8,0.08,'Q_c(1)','fontsize',18)

text(0.8,0.5,'Q_d(1)','fontsize',18)
figure(5)
plot(onebyrclist1,Cc1,onebyrclist1,Cd1,onebyrclist1,Cc0)
xlabel('1/\rho_c')
ylabel('Relative Concentration')
text(0.2,2.4,'C_d(1)','fontsize',18)
text(0.54,1,'C_c(1)','fontsize',18)
text(1,0.5,'C_c(0)','fontsize',18)

function out = des(U,p)
 
y = linspace(0,1,p.N);
dy = y(2)-y(1);

Sd = U(1:p.N);
Qc = U(p.N+1:2*p.N);
Qa = U(2*p.N+1);
Qc1 = U(2*p.N+2);

V = get_conc(U,p);
 
Ss = V(1:p.N);
Qs = V(p.N+1:2*p.N);
Qd = V(2*p.N+1:3*p.N);
Cd = V(3*p.N+1:4*p.N);
Cs= V(4*p.N+1:5*p.N);
Ca =  V(5*p.N+1:6*p.N);
Cc=V(6*p.N+1:7*p.N); 
Fd = Ss./Qs-Sd./Qd;

Fc = -p.DPc+(Sd(p.N)-p.P)./Qc-Ss./Qs;
 


eqSd =[Sd(1)-1,(Sd(2:p.N-1)-Sd(1:p.N-2))/(dy)-p.Hd*Fd(1:p.N-2), ...
    (Sd(p.N)-Sd( p.N-1))/(dy)-p.Hd*Fd(p.N-1)];  
% 
eqQc = [Qc(1)+Qa,(Qc(2:p.N-1)-Qc(1:p.N-2))/(dy)-Fc(1:p.N-2)*p.onebyrc, ...
    (Qc(p.N)-Qc( p.N-1))/(dy)-Fc(p.N-1)*p.onebyrc];

eqQa = Sd(p.N)-1-(Qa+1)*p.rd*p.Hd+p.DPd*p.Hd;
eqQc1 = Qc1-Qc(p.N);

 
out = [eqSd,eqQc,eqQa,eqQc1];

end

function out = get_conc(U,p)

y = linspace(0,1,p.N);
 
Sd = U(1:p.N);
Qc = U(p.N+1:2*p.N);
Qa = U(2*p.N+1);
Qc1 = U(2*p.N+2);

Ss = p.P*(y-1) + Sd(p.N)-Sd;
tm = (1-Sd-p.DPd*p.Hd*y)/(p.rd*p.Hd);
Qs = -1-Qa-Qc+Qc(p.N) - tm ;
Qd = 1+ tm;

Cd=Sd./Qd;
Cs = (p.P+p.DPd*p.Hd).*(1-y)./(Qd+Qa)-p.rd*p.Hd;
Ca=Cd(p.N)-p.P*(y-1)/Qa;
 Cc=-Qa*Ca(1)./Qc;

 out = [Ss,Qs,Qd,Cd,Cs,Ca,Cc];
end
