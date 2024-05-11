
%-------------------------------------------------------------------
%
% Nonlinear compression in the chinchilla cochlea.
% 
% For Chapter 20, Section 20.4.2 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
% Written by James Keener and James Sneyd
%
%-------------------------------------------------------------------

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);

w=linspace(0.6,1.2,500);

mu=0;
w0=1;

a=0.1;
U1=mu^2 + (w0-w).^2;
U2=-mu*mu + 3*(w0-w).^2;
D=27*a*a + 16*(mu^3) - 18*mu*U1;
s=D+(D.*D+4*(U2.^3)).^0.5;
s=s./2;
bigA1= s.^(1/3) + 2*mu - U2./(s.^(1/3));
bigA1 = (bigA1./3).^0.5;

a=0.01;
U1=mu^2 + (w0-w).^2;
U2=-mu*mu + 3*(w0-w).^2;
D=27*a*a + 16*(mu^3) - 18*mu*U1;
s=D+(D.*D+4*(U2.^3)).^0.5;
s=s./2;
bigA2= s.^(1/3) + 2*mu - U2./(s.^(1/3));
bigA2 = (bigA2./3).^0.5;

a=0.0001;
U1=mu^2 + (w0-w).^2;
U2=-mu*mu + 3*(w0-w).^2;
D=27*a*a + 16*(mu^3) - 18*mu*U1;
s=D+(D.*D+4*(U2.^3)).^0.5;
s=s./2;
bigA3= s.^(1/3) + 2*mu - U2./(s.^(1/3));
bigA3 = (bigA3./3).^0.5;

semilogy(w,bigA1)
hold on
semilogy(w,bigA2)
semilogy(w,bigA3)
hold off
xlabel('\omega/\omega_0')
ylabel('A')

%dum=[w' bigA1' bigA2' bigA3'];   % for external plotting
%save('chinchilla.dat','-ascii','dum')

