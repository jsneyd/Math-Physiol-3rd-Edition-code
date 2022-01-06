clear all
close all
clc

%% Response to a bar of width 2R, by direct calculation
a = 2;
R = 5;
lh = 3;

x1 = linspace(-3*R,-R,500);
G1 = (1/(a*a*lh*lh))*abs(exp(-a*abs(x1+R)) - exp(-a*abs(x1-R)));
plot(x1,G1,'LineWidth',2)
hold on

x2 = linspace(-R,R,500);
G2 = (1/(a*a*lh*lh))*( 2 - (exp(-a*abs(x2-R)) + exp(-a*abs(x2+R)) ) );
plot(x2,G2,'LineWidth',2)

x3 = linspace(R,3*R,500);
G3 = (1/(a*a*lh*lh))*abs(exp(-a*abs(x3+R)) - exp(-a*abs(x3-R)));
plot(x3,G3,'LineWidth',2)


%% Response to a bar of width 2R, by convolution with the Green's function

p = @(x) heaviside(x+R).*heaviside(R-x);    % The stimulus. Redefine as needed
gg = @(x) (1/(a*lh*lh))*exp(-a*abs(x));   % The Green's function

x = linspace(-3*R,3*R,500);
dx = x(2)-x(1);

cc = conv(p(x),gg(x),'same')*dx;  % This uses the built-in Matlab convolution function. It's not obvious how to use it to do integrals.

plot(x,cc,'k--','LineWidth',2)