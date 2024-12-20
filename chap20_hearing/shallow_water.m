
%-------------------------------------------------------------------
%
% Simulate shallow-water waves in the cochlea.
% 
% For Chapter 20, Section 20.2.4 of
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
lam = 1.5;
wlist = [800,1500];
for wj=1:2
    w = wlist(wj);     
    L = 3.5;
    l = 0.0035;
    k0 = 1e7;
    r0 = 3000;    
    alpha = sqrt(2*w*w/(l*(k0+1i*w*r0)));
    ar = real(alpha);
    
    % take the square root with positive imaginary part. Be a
    % little careful, as Matlab here gives you a negative imaginary part, which
    % is not the root you want.
    ai = -imag(alpha); 
    
    x = linspace(0,L,2000);
    h = exp(3*lam*x/4 - 2*ai*exp(lam*x/2)/lam + 2*1i*ar*exp(lam*x/2)/lam);
    
    figure(wj)
        plot(x,real(h)/max(abs(h)),'r',x,abs(h)/max(abs(h)),'--b',x,-abs(h)/max(abs(h)),'--b')
        xp = -2*log(4*ai/(3*lam))/lam;
        xlim([0,L])
        xlabel('x (cm)')
        box off
        formatSpecF = '%6.0f\n';
        title(strcat('\omega = ',sprintf(formatSpecF,w),'/s'))
        ylabel('normalized amplitude, Re(h)/|h|_{max}')
end
%igorout = [x' real(eta/max(eta))' abs(eta/max(eta))' -abs(eta/max(eta))'];
%writematrix(igorout,'shallow_800.dat')



