
%  -------------------------------------------------------------------
%
%   Compute solutions of the T-cell model.
%
%   For Chapter 13, Section 13.3 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function T_cell_clone

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global threshold
threshold = false;  % do you want to threshold the bacteria at B<1?  If so, set threshold = true

init = [1 100 0 0];
tspan = linspace(0,40,1000);

warning("off")
[T,Y] = ode15s(@(t,x)rhs(t,x),tspan,init);
B = Y(:,1); N = Y(:,2); A = Y(:,3); M = Y(:,4);
figure(2)

subplot(2,2,1)
plot(T,B,'r')
xlabel('t')
ylabel('B')
subplot(2,2,2)
plot(T,N,'g')
ylabel('N')
xlabel('t')
subplot(2,2,3)
plot(T,A,'b')
xlabel('t')
ylabel('A')
subplot(2,2,4)
plot(T,M,'k')
xlabel('t')
ylabel('M')


figure(1)
plot(T,B)
ylabel('B (number of cells)')
yyaxis right
plot(T,A,T,10*M)
legend('boxoff')
legend('B','A','10 x M')
  
    ylabel('A, 10 x B(numberof cells)')

    xlabel('t (days)')
    
%writematrix([T B N A M],'T_cell_out.dat')   % for external plotting

end % of main

%%
function out=rhs(t,x)
global threshold

if (threshold==true & x(1)<0.9)
    x(1)=0;
end

B = x(1);
N = x(2);
A = x(3);
M = x(4);

sigma = 0;
rN = 0;
dN = 0.001;
aN = 1;
dA = 1;
m = 0.005;
rM = 0;
aM = 0;
dM = 0;

p = 2;
r = 5;
k = 0.001;
h = 1000;

f = B/(h+B);
out(1) = r*B - k*B*A;
out(2) = sigma - dN*N - aN*f*N;
out(3) = f*(aN*N + aM*M + p*A) - dA*A - m*(1-f)*A;
out(4) = m*(1-f)*A - aM*f*M - dM*M;
out = out';

end

