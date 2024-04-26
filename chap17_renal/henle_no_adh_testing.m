%    -------------------------------------------------------------------
%
%     Solve the nephron equations for the loop of Henle (no ADH)
%
%     For Chapter 17, Section 17.3.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function henle_no_adh
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
p.N=151; % number of grid points
p.y = linspace(0,1,p.N);

% make up some initial guess
Qa = -.9;
Qd = -Qa*p.y +(1-p.y);
Sd=-(Qd -1)*(p.rd*p.Hd)+1-p.DPd*p.Hd*p.y ;
X0=[Sd, Qa];
options = optimset('Display','off');
[U]=fsolve(@(x)des(x,p),X0,options);
Qa = U(end)
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




end % of main

%%
function out = des(U,p)
    dy = p.y(2)-p.y(1);
    Sd = U(1:p.N);
    Qa = U(p.N+1);
    [Qd,Qs,Sd,Sa,Ss,Cd,Ca,Cs]=get_fs(U ,p);

    Fd = Ss./Qs-Sd./Qd;

    eqSd =[Sd(1)-1,...
           (Sd(2:p.N-1)-Sd(1:p.N-2))/(dy)-p.Hd*Fd(1:p.N-2), ...
           (Sd(p.N)-Sd( p.N-1))/(dy)-p.Hd*Fd(p.N-1)];
    eqQa = Sd(p.N)-1-(Qa+1)*p.rd*p.Hd+p.DPd*p.Hd;
    out = [eqSd, eqQa ];
    %max(abs(out))
end

%%
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

