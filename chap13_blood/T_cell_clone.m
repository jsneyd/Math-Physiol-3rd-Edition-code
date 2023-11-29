clear all
close all
clc
global threshold
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);
threshold = false;  % do you want to threshold the bacteria at B<1?

init = [1 100 0 0];
tspan = linspace(0,40,1000);

[T,Y] = ode15s(@(t,x)rhs(t,x),tspan,init);
B = Y(:,1); N = Y(:,2); A = Y(:,3); M = Y(:,4);
subplot(2,2,1)
plot(T,B,'r')
subplot(2,2,2)
plot(T,N,'g')
subplot(2,2,3)
plot(T,A,'b')
subplot(2,2,4)
plot(T,M,'k')

writematrix([T B N A M],'T_cell_out.dat')


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

