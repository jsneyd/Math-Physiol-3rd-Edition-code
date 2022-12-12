% This code calculates the fundamental solution and a bunch of Green's
% functions for the cable equation on a finite domain. The Green's functions are calculated
% for sealed-end and closed-circuit ends (i.e., Neumann and Dirichlet BCs)
% and are calculated by a sum of fundamental solutions and by Fourier
% series. So there's a bunch of different curves.

clear all
close all
clc

L = 1;
xi = 0.7*L;
tau = 0;
[T,X] = meshgrid(tau:0.005:tau+0.2,xi-0.8:0.01:xi+0.8);

%%----------------------------------------

% fundamental solution on infinite domain
Z = fund(T,X,xi,tau);
figure(1)
surf(T,X,Z)
xlabel('T')
ylabel('X')
zlabel('V')

%%----------------------------------------

% Green's function on finite domain. Sealed ends.
[T,X] = meshgrid(tau:0.001:tau+0.2,0:L/200:L);
Z1 = GMISE(T,X,xi,tau,L);  % method of images
Z2 = GFSSE(T,X,xi,tau,L);  % Fourier series
figure(2)
subplot(2,1,1)
surf(T,X,Z1)
xlabel('T'); ylabel('X'); zlabel('V')
title('method of images')
subplot(2,1,2)
surf(T,X,Z2)
xlabel('T'); ylabel('X'); zlabel('V')
title('Fourier expansion')
zlim([0,10])

% some plots for fixed T
tf = [0.001 0.01 0.05 0.1];
xf = linspace(0,L,500);
Vf1 = GMISE(tf(1),xf,xi,tau,L);
Vf2 = GMISE(tf(2),xf,xi,tau,L);
Vf3 = GMISE(tf(3),xf,xi,tau,L);
Vf1F = GFSSE(tf(1),xf,xi,tau,L);
Vf2F = GFSSE(tf(2),xf,xi,tau,L);
Vf3F = GFSSE(tf(3),xf,xi,tau,L);
figure(3)
plot(xf,Vf1,'k',xf,Vf2,'k',xf,Vf3,'k','LineWidth',2)
hold on
plot(xf,Vf1F,'r--',xf,Vf2F,'b--',xf,Vf3F,'g--','LineWidth',2)
xlabel('T'); ylabel('V')
legend('method of images','','','Fourier expansion:t=0.001','t=0.01','t=0.05','Location','northwest')
set(gca,'FontSize',14)


%%----------------------------------------

% Green's function on finite domain. Short-Circuit ends.

[T,X] = meshgrid(tau:0.005:tau+0.2,0:L/200:L);
Z1 = GMISC(T,X,xi,tau,L);  % method of images
Z2 = GFSSC(T,X,xi,tau,L);  % Fourier series
figure(4)
subplot(2,1,1)
surf(T,X,Z1)
xlabel('T'); ylabel('X'); zlabel('V')
title('method of images')
subplot(2,1,2)
surf(T,X,Z2)
xlabel('T'); ylabel('X'); zlabel('V')
title('Fourier expansion')
zlim([0,10])

% some plots for fixed T
tf = [0.001 0.01 0.05 0.1];
xf = linspace(0,L,500);
Vf1 = GMISC(tf(1),xf,xi,tau,L);
Vf2 = GMISC(tf(2),xf,xi,tau,L);
Vf3 = GMISC(tf(3),xf,xi,tau,L);
Vf1F = GFSSC(tf(1),xf,xi,tau,L);
Vf2F = GFSSC(tf(2),xf,xi,tau,L);
Vf3F = GFSSC(tf(3),xf,xi,tau,L);
figure(5)
plot(xf,Vf1,'k',xf,Vf2,'k',xf,Vf3,'k','LineWidth',2)
hold on
plot(xf,Vf1F,'r--',xf,Vf2F,'b--',xf,Vf3F,'g--','LineWidth',2)
xlabel('T'); ylabel('V')
legend('method of images','','','Fourier expansion:t=0.001','t=0.01','t=0.05','Location','northwest')
set(gca,'FontSize',14)

%%------------------------------------------
%%----------------------------------------

% Fundamental solution
function out = fund(t,x,xi,tau)
    d1 = 1./sqrt(4*pi*(t-tau));
    d2 = exp(-(x-xi).^2./(4*(t-tau)));
    d3 = exp(-(t-tau));
    out = d1.*d2.*d3;
end

% Base term in Fourier expansion. Sealed ends.
function out = baseF(t,x,xi,tau,n,L)
    d1 = cos(n*pi*xi/L)*cos(n*pi*x/L);
    d2 = exp(-(1+(n*pi/L)^2)*(t-tau));
    out = d1.*d2/L;
end

% Base term in Fourier expansion. Short-circuit ends.
function out = baseFSC(t,x,xi,tau,n,L)
    d1 = sin(n*pi*xi/L)*sin(n*pi*x/L);
    d2 = exp(-(1+(n*pi/L)^2)*(t-tau));
    out = d1.*d2/L;
end

% Greens function, Method of Images, Sealed Ends
function out = GMISE(T,X,xi,tau,L)
    out = fund(T,X,xi,tau) + ...
        fund(T,X,-xi,tau) + fund(T,X,2*L-xi,tau) + ...
        fund(T,X,2*L+xi,tau) + fund(T,X,-2*L+xi,tau);   % Only the first four corrections.
end

% Greens function, Method of Images, Short-Circuit Ends
function out = GMISC(T,X,xi,tau,L)
    out = fund(T,X,xi,tau) - fund(T,X,-xi,tau) + fund(T,X,2*L+xi,tau) - fund(T,X,-2*L-xi,tau) + ...
     fund(T,X,3*L+xi,tau) - fund(T,X,-3*L-xi,tau) ...
     - fund(T,X,2*L-xi,tau) + fund(T,X,-2*L+xi,tau) - fund(T,X,3*L-xi,tau) + fund(T,X,-3*L+xi,tau) - ...
     fund(T,X,4*L-xi,tau);
end

% Greens function, Fourier series, Sealed Ends
function out = GFSSE(t,x,xi,tau,L)
out = baseF(t,x,xi,tau,0,L);
for i=1:20
    out = out + baseF(t,x,xi,tau,i,L) + baseF(t,x,xi,tau,-i,L);
end
end

% Greens function, Fourier series, short-ciruit Ends
function out = GFSSC(t,x,xi,tau,L)
out = baseFSC(t,x,xi,tau,0,L);
for i=1:20
    out = out + baseFSC(t,x,xi,tau,i,L) + baseFSC(t,x,xi,tau,-i,L);
end
end

