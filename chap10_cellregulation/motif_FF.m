clear all
close all
clc

init=[0 0];
tspan=linspace(0,50,5000);

[T,Y] = ode15s(@rhs,tspan,init);
R = Y(:,1);
X = Y(:,2);
S = heaviside(T-10) + heaviside(T-20) + heaviside(T-30);

plot(T,R,T,S)

out = [T R S];
save('test.dat','out','-ascii')

%%
function out = rhs(t,Y)
R = Y(1);
X = Y(2);

n=2;
K=2;
S = heaviside(t-10) + heaviside(t-20) + heaviside(t-30);
phi = K*X^n/(K+X^n);
out(1) = S - phi*R;
out(2) = S - X;
out = out';

end