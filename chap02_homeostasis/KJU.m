%-------------------------------------------------------------------

% Matlab code for solving the KJU model of a Na-transporting cell. The
% nonlinear equations are solved using a built-in Matlab nonlinear solver, fsolve.

% For Chapter 2, Fig. 2.21 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------
function KJU
clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);
Nm = linspace(10,210,50);
x0 = [20 100 0.1 -2];
 RTbyF = 25.8;
for i=1:50
x = fsolve(@(x)KJUfun(x,Nm(i)),x0);
x0 = x;
Ni(i) = x(1);
Ki(i) = x(2);
mu(i) = x(3);
v(i) = x(4);
end

figure(2)
plot(Nm,mu)
xlabel('mucosal [Na^+], N_m (mM)')
ylabel('\mu (scaled cell volume)')
figure(1)
plot(Nm,RTbyF*v)
xlabel('mucosal [Na^+], N_m (mM)')
ylabel('V (mV)')
ylim([-80,0])
end

%%
function out = KJUfun(x,Nm)
  Ks = 2.5;
  Cs = 122.5;
  Ns = Cs - Ks;
  z = -2;
  rhon = 1.5;
  rhok = 1;
  
  Ni = x(1);
  Ki = x(2);
  mu = x(3);
  v = x(4);

  dum = exp(-v);
  Ci = Cs/dum;
    
  out(1) = v*(Ni-Nm*dum)/(1-dum) + 3*rhon*Ni;
  out(2) = v*(Ki-Ks*dum)/(1-dum) - 2*rhok*Ni;
  out(3) = (Ni+Ki-Ci) + z/mu;
  out(4) = Ni+Ki+Ci + 1/mu - (Ns+Ks+Cs);
end