
% Hai-Murphy-Huxley model of smooth muscle. Compute the maximal velocity
% and the isometric force

function muscle_smooth_max_v
close all
clear all
clc

par.f1 = 0.88; 
par.g1 = 0.21;
par.g2 = 4.4;
par.gL1 = 0.01;
par.gL2 = 0.2;
par.h=1;
par.k1 = 0.35;
par.k2 = 0.1;
par.k5 = 0.35;
par.k6 = 0.1;

par.num = 250;      % number of space points
par.numt = 150;     % number of time outputs
par.tend = 2;       % final time
par.v = 0.5;          % velocity

get_isometric_force(par);
vmax = fzero(@(v)get_load(v,par),2);
fprintf('maximal velocity is %5.4f\n',vmax) % The maximal velocity occurs when the load = 0

% now calculate and plot the distributions at the maximal velocity
par.v = vmax;
xspan = linspace(4,-6,2000);
n0 = [1,0,0,0];  % order is Nm Nmp Nam Namp
[xout,nout]=ode15s(@(x,n)derivs_steady(x,n,par),xspan,n0);
Nm = nout(:,1);
Nmp = nout(:,2);
Nam = nout(:,3);
Namp = nout(:,4);

plot(xout,Nam,'b','LineWidth',2)
hold on
plot(xout,Namp,'r','LineWidth',2)
xlabel('x')
legend('Nam','Namp')

end

%%
function load = get_load(v,par) % Compute the force for a given velocity
xspan = linspace(2,-6,2000);
n0 = [1,0,0,0];  % order is Nm Nmp Nam Namp
par.v = v;
[xout,nout]=ode15s(@(x,n)derivs_steady(x,n,par),xspan,n0);
Nam = nout(:,3);
Namp = nout(:,4);
load = -trapz(xout,xout.*(Nam + Namp));  % Has to be minus since x is going backwards in the integration
end

%%
function get_isometric_force(par) % Compute the isometric steady-state force. Just for fun.
x0 = linspace(-10,10,par.num);
k1=par.k1;
k2=par.k2;
k5=par.k5;
k6=par.k6;
k7=gL(x0,par);
steady = (f(x0,par)*k1.*(k5 + k7))./(f(x0,par)*k1*k5 + f(x0,par)*k1*k6 + f(x0,par)*k1.*k7 + f(x0,par)*k6.*k7 + g(x0,par)*k1*k5 + g(x0,par)*k2*k5 + g(x0,par)*k1.*k7 + g(x0,par)*k2.*k7 + k1*k6*k7 + k2*k6*k7) ...
            + (f(x0,par)*k1*k6)./(f(x0,par)*k1*k5 + f(x0,par)*k1*k6 + f(x0,par)*k1.*k7 + f(x0,par)*k6.*k7 + g(x0,par)*k1*k5 + g(x0,par)*k2*k5 + g(x0,par)*k1.*k7 + g(x0,par)*k2.*k7 + k1*k6*k7 + k2*k6*k7);
fprintf('isometric force is %5.4f\n',trapz(x0,x0.*steady))
end

%% RHS of ode for computing the steady-state distributions
function out = derivs_steady(x,n,par)
% Solve these odes backwards in x to get the steady-state distributions. 
nm = n(1);
nmp = n(2);
nam = n(3);
namp = n(4);

out(1) = -(1/par.v)*(par.k2*nmp - par.k1*nm + gL(x,par).*nam);                       % Nm
out(2) = -(1/par.v)*(par.k1*nm - (par.k2 + f(x,par)).*nmp + g(x,par).*namp);     % Nmp
out(3) = -(1/par.v)*(par.k6*namp - (par.k5+gL(x,par)).*nam);                         % Nam
out(4) = -(1/par.v)*(par.k5*nam + f(x,par).*nmp - (par.k6+g(x,par)).*namp);      % Namp
out = out';
end

%% ------------
function out = f(x,par)
out = 0 + (x>0 & x<=1).*par.f1.*x.*(1-x);
end

function out = g(x,par)
out = (x<0).*par.g2 + (x>0).*par.g1.*x;
end

function out = gL(x,par)
out = (x<0).*par.gL2 + (x>0).*par.gL1.*x;
end





