
%  -------------------------------------------------------------------
%
%   ODE integrator for the YNI equations.
%
%   For Chapter 12, Section 12.2.1 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  The vector S = [ V m h j d f x Cai] contains                             %
%  V (the membrane potential), and gating variables f, h, p, q, d.                                                        %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function YNI_ode

close all
clear all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);
 %The natural period is about 800'
 global gnabar Vna
      gnabar = 0.5;
      Vna = 30.;

Inp = 0;
BCL = 25;

tstep = 1;
t_end = 1000;

V = -43.2572;
%V = -60;

f = 1;
h = 0.9207;
p = .1003;
q = 0.0365;
d = 0.0001;

%specify the output points
tspan = [0:tstep:t_end];

s0 = [V;f;h;p;q;d];

[T,S] = ode15s(@deRHS,tspan, s0, odeset('maxstep',1));
V = S(:,1);
f = S(:,2);
h = S(:,3);
p = S(:,4);
q = S(:,5);
d = S(:,6);


figure(1)
 
plot(T,V,'linewidth',2)
 xlabel('t (ms)')
 ylabel('V (mV)')

 box off
 figure(2)
plot(T,f,'b',T,h,'r',T,p,'g',T,q,T,d,'linewidth',2)
legend('f','h','p','q','d','fontsize',16)
 xlabel('t (ms)')

 for j = 1:length(T)
      currents= IV(S(j,:));
      INa(j) = currents(1);
      IK(j) = currents(2);
      Ix(j) = currents(3);
      Is(j) = currents(4);
      Il(j) = currents(5);
 end

figure(3);
 plot(T,INa,'k',T,IK,'r',T,Ix,'b',T,Is,'m',T,Il,'linewidth',2);
xlabel('t (ms)');

legend('I_{Na}','I_K','I_h','I_s','I_l', 'fontsize',16)

end % of main


%%
%the right hand side for ode simulation:
function s_prime=deRHS(t,s)
    global gnabar Vna
    V = s(1);
    f = s(2);
    h = s(3);
    p = s(4);
    q = s(5);
    d = s(6);

    %Remark:  I think there is a mistake in ad

    %Remark: q is the slowest variable.
    %YNI model:
    %'The natural period is about 300'

    ah = 0.001209*exp(-(V+20.)/6.534);
    bh = 1./(exp(-0.1*(V +30.)) + 1.);
    ap = 0.009./(1.+exp(-(V+3.8)/9.71)) + 0.0006;
    bp = 0.000225*(V+40.)./(exp((V+40.)/13.3) -1.);
    aq = 0.00034*(V+100.)./(exp((V+100.)/4.4)-1.) +.0000495;
    bq = 0.0005*(V+40.)./(1-exp(-V/6-6.6667)) + .0000845;
    tm = V+35.;
    ad=(tm>0).*tm./(1.-exp(-tm/2.5))+...
    (tm<0).*tm.*exp(tm/2.5)./(exp(tm/2.5)-1.)+...
    (tm==0).* 2.5;
    ad = (V>0).*(.01045*ad + .03125*V./(1-exp(-V/4.8)))+...
    (V<0).* 0.01045.*ad + .03125*V.*exp(V/4.8)./(exp(V/4.8)-1.)+...
    (V==0).*( .01045*ad +.03125*4.8);

    bd = .00421*(V-5.)./(exp(V/2.5-2.)-1.);
    af = .000355*(V+20.)./(exp((V+20.)/5.633)-1.);
    bf = .000944*(V+60.)./(1.+exp(-(V+29.5)/4.16));

    Fs = af*(1.-f) - bf*f;
    Fh = ah*(1.-h) - bh*h;
    Fp = ap*(1.-p) - bp*p;
    Fq = aq*(1.-q) - bq*q;
    Fd = ad*(1.-d) - bd*d;

    currents = IV(s);
    INa= currents(1);
    IK = currents(2);
    Ih = currents(3);
    Is = currents(4);
    Il = currents(5);

    Fv = -(Is + INa + IK + Ih + Il);
    %Fv = 0;
    s_prime = [Fv Fs Fh Fp Fq Fd]';

end

%%
function currents=IV(s)
    global gnabar Vna
    V = s(1);
    f = s(2);
    h = s(3);
    p = s(4);
    q = s(5);
    d = s(6);

    tm = V+37.;
    am = (tm>0).*tm./(1.-exp(-tm/10.))+...
    (tm<0).*tm.*exp(tm/10.)./(exp(tm/10.)-1.)+...
    (tm==0)*10;
    bm = 40.*exp(-0.056*(V+62.));
    minf = am./(am+bm);
    INa = gnabar*minf.^3.*h.*(V-Vna);
    Ik = 0.7*p*(exp(0.0277*(V+90.))-1.)/(exp(0.0277*(V+40.)));
    Ih = 0.4*q*(V+25.);
    Is = 12.5*(0.95*f+0.05)*(exp(V/15.-2.)-1.)*(0.95*d+0.05);
    Il = 0.8*(1.-exp(-V/20.-3.));

    currents = [INa, Ik, Ih, Is, Il];
end
