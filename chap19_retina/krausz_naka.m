% Solutions of the Krausz-Naka model of photoreceptor cell/horizontal cell
% interactions in the catfish retina. This only computes responses to
% steady inputs.

clear all
close all
clc

%% Response to a bar of width 2R, by direct calculation

% Usually alpha is a function of w, but here we're just setting it to be constant, to determine the resopnse to steady inputs.
% Just because it's easier.

alpha = 3.77;  
R = 2;
lh = 0.267;

x1 = linspace(-3*R,-R,500);
G1 = (1/(alpha*alpha*lh*lh))*abs(exp(-alpha*abs(x1+R)) - exp(-alpha*abs(x1-R)));
plot(x1,G1,'LineWidth',2)
hold on

x2 = linspace(-R,R,500);
G2 = (1/(alpha*alpha*lh*lh))*( 2 - (exp(-alpha*abs(x2-R)) + exp(-alpha*abs(x2+R)) ) );
plot(x2,G2,'LineWidth',2)

x3 = linspace(R,3*R,500);
G3 = (1/(alpha*alpha*lh*lh))*abs(exp(-alpha*abs(x3+R)) - exp(-alpha*abs(x3-R)));
plot(x3,G3,'LineWidth',2)


%% Response to a bar of width 2R, by convolution with the Green's function

p = @(x) heaviside(x+R).*heaviside(R-x);    % The stimulus. Redefine as needed
gg = @(x) (1/(alpha*lh*lh))*exp(-alpha*abs(x));   % The Green's function

x = linspace(-3*R,3*R,500);
dx = x(2)-x(1);

cc = conv(p(x),gg(x),'same')*dx;  % This uses the built-in Matlab convolution function. It's not obvious how to use it to do integrals.

plot(x,cc,'k--','LineWidth',2)