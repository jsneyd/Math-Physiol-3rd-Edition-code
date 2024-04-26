%    -------------------------------------------------------------------
%
%     Solve the nephron equations for the loop of Henle.
%
%     For Chapter 17, Section 17.3.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function henle
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

onebyrclist0 = [0.,0.5];
options = optimset('Display','off');

for j = 1:2
    p.onebyrc=onebyrclist0(j);

    % make up some initial guess

    % The problem is extremely sensitive to the initial condition, and
    % considerable care is needed, particularly in the choice of Qa. It is
    % easy to find an initial value for Qa that works in Matlab but
    % doesn't work in Octave.

    Qa = -.22;
    Qd = -Qa*y +(1-y);
    Sd=-(Qd -1)*(p.rd*p.Hd)+1-p.DPd*p.Hd*y ;
    Qc1=0.1;
    Qc = Qc1*y -Qa*(1-y);

    X0=[Sd,Qc,Qa,Qc1];

    [U]=fsolve(@(x)des(x,p),X0,options);
    Sd = U(1:p.N);
    Qc = U(p.N+1:2*p.N);
    Qa = U(2*p.N+1);
    Qc1 = U(2*p.N+2);
    [Ss,Qs,Qd,Cd,Cs,Ca,Cc] = get_conc(U,p);

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

    if (j==1)
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
            plot(y,Cd,y,Ca,'--',y,Cc)
            legend('C_d','C_a','C_c','location','northwest')
            xlabel('y')
            ylabel('Relative Concentration')
            box off
    end
 end

% now loop on rc for no ADH case:
% make up some initial guess
Qa = -.5;
Qd = -Qa*y +(1-y);
Sd=-(Qd -1)*(p.rd*p.Hd)+1-p.DPd*p.Hd*y ;
Qc1=0.1;
Qc = Qc1*y -Qa*(1-y);
X0=[Sd,Qc,Qa,Qc1];
onebyrclist1 = [0.01:0.025:2];

for k = 1:length(onebyrclist1)
    p.onebyrc = onebyrclist1(k);
    [U]=fsolve(@(x)des(x,p),X0,options);
    % use the solution as the initial guess for the next try
    X0=U;
    Sd = U(1:p.N);
    Qc = U(p.N+1:2*p.N);
    Qa = U(2*p.N+1);
    Qc1 = U(2*p.N+2);

    [Ss,Qs,Qd,Cd,Cs,Ca,Cc] = get_conc(U,p);

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
    legend('C_c(1)','C_d(1)','C_c(0)')

end % of main

%%
function out = des(U,p)
    y = linspace(0,1,p.N);
    dy = y(2)-y(1);

    Sd = U(1:p.N);
    Qc = U(p.N+1:2*p.N);
    Qa = U(2*p.N+1);
    Qc1 = U(2*p.N+2);

    [Ss,Qs,Qd,Cd,Cs,Ca,Cc] = get_conc(U,p);

    Fd = Ss./Qs-Sd./Qd;
    Fc = -p.DPc+(Sd(p.N)-p.P)./Qc-Ss./Qs;

    eqSd =[Sd(1)-1,...
           (Sd(2:p.N-1)-Sd(1:p.N-2))/(dy)-p.Hd*Fd(1:p.N-2), ...
           (Sd(p.N)-Sd( p.N-1))/(dy)-p.Hd*Fd(p.N-1)];
    %
    eqQc = [Qc(1)+Qa,...
            (Qc(2:p.N-1)-Qc(1:p.N-2))/(dy)-Fc(1:p.N-2)*p.onebyrc, ...
            (Qc(p.N)-Qc( p.N-1))/(dy)-Fc(p.N-1)*p.onebyrc];
    eqQa = Sd(p.N)-1-(Qa+1)*p.rd*p.Hd+p.DPd*p.Hd;
    eqQc1 = Qc1-Qc(p.N);

    out = [eqSd,eqQc,eqQa,eqQc1];
end

%%
function [Ss,Qs,Qd,Cd,Cs,Ca,Cc] = get_conc(U,p)
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
    % using a na-dependent ATPase:
    %Ca = Cd(p.N)*sqrt(2*Qa./(p.P*Cd(p.N)^2.*(y-1)+2*Qa));
    Cc=-Qa*Ca(1)./Qc;
end
