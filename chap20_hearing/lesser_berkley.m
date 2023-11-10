clear all
close all
clc

num = 1000;
L = 3.5;   % Units of cm
l = 0.35;
lam = 1.5;
x = linspace(0,1,num);
k = (1e7)*exp(-lam*x*L);   % don't forget that x is scaled by L
mass = 0.05;
r = 3000*exp(-lam*x*L);
w = 800;
Z = 1i*w*mass + r + k/(1i*w);
W = Z/(1i*w*L);
sigma = l/L;


% Initialise
N =150;
alpha = zeros(N+1);
f = zeros(N+1,1);
A = zeros(N+1,1);

% get matrix coefficients
for n=0:N
    nn = n+1;
    for m = 0:N
        mm = m+1;
        alpha(mm,nn) = 2*cosh(n*pi*sigma)*trapz(x,cos(n*pi*x).*cos(m*pi*x)./W);
        f(mm) = -trapz(x,x.*(2-x).*cos(m*pi*x)./W);
    end
    alpha(nn,nn) = alpha(nn,nn) + n*pi*sinh(n*pi*sigma)/2;
end
f(1) = f(1) - sigma;

% solve linear equations for Fourier coefficients
A  = alpha\f;

% construct phi on y=0
y = 0.0;   % evaluate eta on the membrane
phi = x.*(1-x/2) - sigma*y.*(1-y/(2*sigma));
for n=0:N
    phi = phi + A(n+1)*cosh(n*pi*(sigma-y))*cos(n*pi*x);
end

Fhat = 1;   % driving force
phi = phi/(1i*w*L*Fhat);
eta = 2*phi./W;

figure(2)
plot(x,eta,'r','LineWidth',1)
hold on
plot(x,abs(eta),'--b',x,-abs(eta)','--b','LineWidth',2)


% animate the wave, just for funsies
N = 200;
times = linspace(0,10/w,N);
figure(3)
for i=1:N
    etawave = eta*exp(1i*w*times(i));
    plot(x,etawave,'r',x,abs(etawave),'--b',x,-abs(etawave),'--b')
    drawnow
    pause(0.02)
end


%igorout = [x' real(eta/max(eta))' abs(eta/max(eta))' -abs(eta/max(eta))'];
%writematrix(igorout,'igor_800.dat')



