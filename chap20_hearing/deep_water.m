%-------------------------------------------------------------------
%
% Simulate deep-water waves in the cochlea.
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
L = 3.5;   % Units of cm
l = 3.5;
N = 1000;
x = linspace(0,L,N);
wlist = [800,1500];
igorout = [];
for wj=1:2
    w=wlist(wj);    
    sigma = l/L;    
    k0 = 1e7;
    r0 = 3000;
    F0 = 1000;
    lam = 1.5;
    rho = 1;
    k = k0*exp(-lam*x); 
    r = r0*exp(-lam*x);
    Z = (r + k/(1i*w));
    Z0 = (r0 + k0/(1i*w));
    Y = 2./Z;
    h0 = F0/Z0;
    
    intY = 2*1i*w*(exp(lam*x) - 1)/(lam*(1i*w*r0 + k0));
    h = -1i*rho*h0.*Y.*exp(-w*rho*intY);
       
    figure(wj)
        plot(x,h/max(abs(h)),'r')
        hold on
        plot(x,abs(h)/max(abs(h)),'--b','LineWidth',2)
        plot(x,-abs(h)/max(abs(h)),'--b','LineWidth',2)
        hold off
        xlim([0,L])
        box off
        formatSpecF = '%6.0f\n'; 
        title(strcat('\omega = ',sprintf(formatSpecF,w),'/s'))
        ylabel('normalized amplitude')
        xlabel('x (cm)')

    igorout = [igorout x' real(h/max(abs(h)))' abs(h/max(h))' -abs(h/max(h))'];
end

%writematrix(igorout,'deep.dat')  % for external plotting




% figure(2)
% % test with the analytic solution
% beta = 2*rho*(w^2)*(w*r0 + 1i*k0)/(lam*(k0^2 + w^2*r0^2));
% xi = w*rho/(1i*w*r0+k0);
% eta2 = 2*eta0*xi*exp(lam*x + beta*(1-exp(lam*x)));
% plot(x,eta2,'r','LineWidth',2)
% hold on
% plot(x,abs(eta2),'b','LineWidth',2)



