clear all
close all
clc

lam = 1.5;
w = 800;
L = 3.5;
l = 0.0035;

k0 = 1e7;
r0 = 3000;

alpha = sqrt(2*w*w/(l*(k0+1i*w*r0)));
ar = real(alpha);

% take the square root with positive imaginary part. You have to be a
% little careful, as Matlab here gives you a negative imaginary part, which
% is not the root you want.
ai = -imag(alpha); 

x = linspace(0,L,2000);
eta = exp(3*lam*x/4 - 2*ai*exp(lam*x/2)/lam + 2*1i*ar*exp(lam*x/2)/lam);

figure(2)
plot(x,eta/max(eta),'r',x,abs(eta)/max(abs(eta)),'--b',x,-abs(eta)/max(abs(eta)),'--b')
xp = -2*log(4*ai/(3*lam))/lam
xlim([0,L])

%igorout = [x' real(eta/max(eta))' abs(eta/max(eta))' -abs(eta/max(eta))'];
%writematrix(igorout,'shallow_800.dat')



