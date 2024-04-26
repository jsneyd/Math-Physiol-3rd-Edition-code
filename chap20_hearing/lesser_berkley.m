% code to simulate cochlear waves with the Lesser-Berkley model
% 
% % For Figure  20.7 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

% Written by James Keener and James Sneyd

%-------------------------------------------------------------------
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);
num = 200;
L = 3.5;   % Units of cm
l = 0.35;
lam = 1.5;
xi = linspace(0,1,num);
k = (1.e7)*exp(-lam*xi*L);   % don't forget that  x=xi*L  
mass = 0.05;
r = 3000*exp(-lam*xi*L);
% two cases:
wlist=[800,1500];

for wj = 1:2
    w = wlist(wj);


Z = 1i*w*mass + r -1i* k/w;
W = -1i*Z/(w*L);
sigma = l/L;


% Initialize
N =100;
alpha = zeros(N+1);
f = zeros(N+1,1);
A = zeros(N+1,1);

% get matrix coefficients
for n=0:N
    nn = n+1;
    for m = 0:N
        mm = m+1;
        alpha(mm,nn) = 2*cosh(n*pi*sigma)*trapz(xi,cos(n*pi*xi).*cos(m*pi*xi)./W);
        f(mm) = -trapz(xi,xi.*(2-xi).*cos(m*pi*xi)./W);
    end
    alpha(nn,nn) = alpha(nn,nn) + n*pi*sinh(n*pi*sigma)/2;
end
f(1) = f(1) - sigma;

% solve linear equations for Fourier coefficients
A  = alpha\f;

 
% construct phi on y=0
y = 0.0;   % evaluate eta on the membrane
phi = xi.*(1-xi/2) - sigma*y.*(1-y/(2*sigma));
for n=0:N
    phi = phi + A(n+1)*cosh(n*pi*(sigma-y))*cos(n*pi*xi);
end

Fhat = 1;   % driving force
phi = phi/(1i*w*L*Fhat);
eta = 2*phi./W;

figure(2*wj-1)
plot(xi*L,real(eta),'r' )
hold on
plot(xi*L,abs(eta),'--b',xi*L,-abs(eta)','--b','LineWidth',2)
 xlabel('x (cm)')
   ylabel('amplitude')
    xlim([0,L])
   formatSpecF = '%6.0f\n';

   title(strcat('\omega = ',sprintf(formatSpecF,w),'/s'))
   box off

% animate the wave, just for fun 
N = 200;
times = linspace(0,10/w,N);
figure(2*wj)
for j=1:N
    etawave = eta*exp(1i*w*times(j));
    plot(xi*L,real(etawave),'r',xi*L,abs(etawave),'--b',xi*L,-abs(etawave),'--b')
   xlabel('x (cm)')
   ylabel('amplitude')
   box off
    formatSpecF = '%6.0f\n';

   title(strcat('\omega = ',sprintf(formatSpecF,w),'/s'))
   xlim([0,L])
   drawnow
    pause(0.02)
end

end


%igorout = [x' real(eta/max(eta))' abs(eta/max(eta))' -abs(eta/max(eta))'];
%writematrix(igorout,'igor_800.dat')



