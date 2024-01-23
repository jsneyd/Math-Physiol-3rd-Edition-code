% this is an attempt to solve the nephron equations

% parameters
p.P = 0.9;
p.DPd = 0.15;
p.DPc = 0.22;
p.Hd = 0.1;
p.rd = 0.15;
p.rc = 2.0;

p.N=50; % number of grid points
% make up some initial guess

X0=[];
fsolve(@(x)des(X,p),X0)

function out = des(U,p)
y = linspace(0,1,p.N)
dy = y(2)-y(1);

Sd = U(1:p.N);
Qc = U(p.N+1,2*p.N);
Qa = U(2*p.N+1);
Qc1 = U(2*p.N+2);

Ss = p.P*(y-1) + Sd(p.N)-Sd;
tm = (1-Sd-p.DPd*Hd*y)/(p.rd*p.Hd);
Qs = -1-Qa-Qc+Qc(p.N) - tm ;
Qd = 1+ tm;
Fd = Ss./Qs-Sd./Qd;
Fc = -DPc+Sc(1)./Qc-Ss./Qs;

Cd=(1+p.rd*p.Hd(1-Qd)-p.DPd*p.Hd*y)./Qd;
Cs = (p.P+p.DPd*p.Hd).*(11-y)./(Qd+Qa)-p.rd*p.Hd;


eqSd =[Sd(1)-1,(Sd(3:N)-Sd(1:p.N-2))/(2*dy)-p.Hd*Fd(2:p.N-1), ...
    (Sd( N)-Sd( p.N-1))/(dy)-p.Hd*Fd(p.N)];  
% 
eqQc = [Qc(1)-1-(1-Sd(1)-p.DPd*p.Hd)/(p.rd*p.Hd),(Qd(3:p.N)-Qd(1:p.N-2))/(2*dy)-Fc(2:N-1)/p.rc, ...
    (Qd(p.N)-Qd( p.N-1))/(dy)-Fc(N )/p.rc];

eqQa = Sd(p.N)-1-(Qa+1)*p.rd*p.Hd+p.Drd*p.Hd;
eqQc1 = Qc1-Qc(p.N);

 
out = [eqSd,eqQc,eqQa,eqQc1];

end

