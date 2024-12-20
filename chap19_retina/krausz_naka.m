
% ---------------------------
% Solutions of the Krausz-Naka model of photoreceptor cell/horizontal cell
% interactions in the catfish retina. 
%
% Used to generate the image in Fig. 17 of Chapter 19 of
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
% Written by James Keener and James Sneyd
% ---------------------------

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);
% Work in time units of milliseconds
lh = 0.3;                      % The space constant of electrical diffusion in the horizontal cell layer
A = 3000;
tau = 25 ; 
t0 = 0.02; 

%% Set up and plot k and khat

% First, set up k and khat. We do this by taking the inverse Fourier
% transform of the function k(t), which is given in the text. We plot k(t)
% just for fun, to make sure it looks OK

Fs = 2;              % sampling frequency
t = 0:1/Fs:3000;             % time vector
L = length(t);              % signal length

k = (3/tau)*exp(-(t-t0)/tau).*(1 - exp(-(t-t0)/tau) ).^2;   % From the Krausz-Naka paper.
figure(1)
    plot(t,k)
    ylabel('k(t)')
    xlabel('t')
box off
n = 2^nextpow2(L);          % Extend sample length and pad with zeros
khat = fft(k,n)/n;          % Take the FFT.
f = Fs*(0:(n/2))/n;         % Define the frequencies for the frequency domain
figure(2)
    plot(f,abs(khat(1:n/2+1)))
    %ylabel('\hat{k}(t)')
    xlabel('\omega (Hz)')
box off
% Once we know the FT of k, we can then calculate alpha(w).
alpha = (1 + A*khat).^0.5/lh;
alpha0 = (1+A*khat(1))^0.5/lh;      % The zero frequency value of alpha
fullfield = 1/(alpha0^2*lh^2)       % The full-field response. Used as a check, mostly.

%% Computation of the frequency response of the bar/field ratio at x=0
R = 0.1;           % Set the width of the bar. R=1 is a wide bar. R=0.1 is a thin bar
x = 0;
F = (2 - exp(-alpha.*(x+R)) - exp(alpha.*(x-R)));  % F as defined in the text
ratio = F/2;

figure(3)
    loglog(f,abs(ratio(1:n/2+1)))
    xlabel('\omega (Hz)')
    ylabel('bar/field response ratio')
    box off
%% Response to a steady bar of width 2R, by direct calculation and then by convolution. Plot on same graph, to compare

% To get a steady response, choose alpha to be the value at 0 frequency. So
% use alpha0 here. 

upper = 2*R;
lower = upper;

% Calculate the bar response in three pieces, To the
% left of the bar, inside the bar, and then to the right of the bar.
x1 = linspace(-lower,-R,500);
G1 = (1/(2*alpha0*alpha0*lh*lh))*abs(exp(-alpha0*abs(x1+R)) - exp(-alpha0*abs(R-x1)));
figure(4)
    plot(x1,G1,'LineWidth',2)
    hold on
    
    x2 = linspace(-R,R,500);
    G2 = (1/(2*alpha0*alpha0*lh*lh))*( 2 - (exp(alpha0*(x2-R)) + exp(-alpha0*(x2+R)) ) );
    plot(x2,G2,'LineWidth',2)
    box off
    x3 = linspace(R,upper,500);
    G3 = (1/(2*alpha0*alpha0*lh*lh))*abs(exp(-alpha0*abs(x3+R)) - exp(-alpha0*abs(R-x3)));
    plot(x3,G3,'LineWidth',2)
    ylabel('Response to a steady bar of width 2R')
    xlabel('x')

% Now find the response by convolution
        
p = @(x) heaviside(x+R).*heaviside(R-x);                % The stimulus. Redefine as needed
gg = @(x) (1/(2*alpha0*lh*lh))*exp(-alpha0*abs(x));     % The Green's function
x = linspace(-2*lower,2*upper,500);                     % Increase the bounds, to try and avoid edge effects
dx=x(2)-x(1);

% Use the built-in Matlab convolution function, which is faster than doing the integral directly.
% But be careful, as this can generate edge effects (as it doesn't always correspond to an integral on
% a domain that is large enough), and it also becomes weird as R gets large; it 
% doesn't pick up sharp edges very well, and the edges get sharper as R
% gets larger.
cc = conv(p(x),gg(x),'same')*dx;  

plot(x,cc,'k--','LineWidth',2)
hold off

%% Check
% Now check that the ratio of the bar/field response (computed directly by
% solving the differential equations) is the same as that given by the
% frequency response

bar_response = G2(250);  % Take the midpoint of the bar response, to get the value at x=0.
fprintf('direct solution gives the bar/field ratio as %d\n\n',bar_response/fullfield);
fprintf('the frequency response gives the bar/field ratio as %d\n\n',abs(ratio(1)));






