
%  -------------------------------------------------------------------
%
%   The model of fetal circulation.
%
%   For Chapter 11, Section 11.7 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd.
% 
%  ------------------------------------------------------------------- 

clear all
close all
clc

syms Psa Psv Ps Rsa Rsv Ppa Ppv 
syms Cld Crd Csv Csa Cp
syms Vsv V0p Vp V0s Vsa Vt
syms Rd Rp F 
syms Qr Qs Qd Qp Ql Qf

% Simplify the equations
Crs = 0;
Cls = 0;

% Define the equations
eq1 = Psa - Ps - Qs * Rsa;
eq2 = Ps - Psv - Qs * Rsv;
eq3 = Ppa - Ppv - Qp * Rp;
eq4 = F * (Cld * Ppv - Cls * Psa) - Ql;
eq5 = F * (Crd * Psv - Crs * Ppa) - Qr;
eq6 = Ql + Qd - Qs;
eq7 = Qd + Qp - Qr;
eq8 = Qr + Qf - Qs;
eq9 = Csv * (Psv + Ps) / 2 - Vsv;
eq10 = Cp * (Ppa + Ppv) / 2 + V0p - Vp;
eq11 = Csa * (Psa + Ps) / 2 + V0s - Vsa;
eq12 = Ppa - Psa - Qd * Rd;
eq13 = Vt - Vsa - Vsv - Vp;

eqns = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13];


%% The foramen-open solution
fprintf('The foramen-open solution\n\n')
eqns_open = subs(eqns,Psv,Ppv);
sol = solve(eqns_open, ...
    [Psa, Ppa, Ppv, Ps, Vp, Vsv, Vsa, Vt, Qf, Qs, Qd, Qp, Ql]);
pretty(simplify(sol.Qs,150))
pretty(simplify(sol.Qp,150))
pretty(simplify(sol.Ql,150))


%% The foramen-closed solution
fprintf('\n\nThe foramen-closed solution\n\n')
eqns_closed = subs(eqns,Qf,0);
sol = solve(eqns_closed, ...
    [Psa, Ppa, Ppv, Psv, Ps, Vp, Vsv, Vsa, Vt, Qs, Qd, Qp, Ql]);
pretty(simplify(sol.Ql,150))
pretty(simplify(sol.Ppv,150))
pretty(simplify(sol.Ppv/sol.Psv,150))


