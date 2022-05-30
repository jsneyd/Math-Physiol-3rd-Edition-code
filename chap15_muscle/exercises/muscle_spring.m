
% Huxley model of muscle pulling against spring. 

% We don't actually need to include an ode for the spring extension (L), but we do anyway (it's the last component of y),
% as this provides a useful check on whether or not your integration is accurate.
% function muscle_char_spring

function muscle_spring
close all
clear all
clc

par.f1=43.3; 
par.g1=10;
par.g2=209;
par.h=1;
par.num = 250; % number of space points
numt = 150;    % number of time outputs
tend = 0.2;   % final time

% initial conditions
x0 = linspace(-2,2,par.num);
n0 = zeros(1,par.num);
L0 = 0;                 % The extension of the spring
y0 = [x0 n0 L0];

output_times = linspace(0,tend,numt);  % to make the curves for a constant v
% integrate the odes, one for each x point
[tout,yout]=ode15s(@(t,y)derivs(t,y,par),output_times,y0);

% All the real work is done now. The rest is just plotting the output.
figure(1)
hold on
% calculate velocity (and plot n) for each time
for i=1:numt
    xx = yout(i,1:par.num);
    nn = yout(i,par.num+1:2*par.num);
    v(i) = trapz(xx, xx.*(f(xx,par).*(1-nn) - g(xx,par).*nn))/(1+trapz(xx,nn));
    plot(xx,nn)
end
hold off

figure(2)
L = yout(:,2*par.num+1);
[hAx,hLine1,hLine2] = plotyy(tout,v,tout,L);  % plot the spring extension and velocity over time
xlabel('time')
ylabel(hAx(1),'v')
ylabel(hAx(2),'L')
set(hLine1,'LineWidth',2);
set(hLine2,'LineWidth',2);

fprintf('initial velocity = %5.4f\n', trapz(x0, x0.*(f(x0,par).*(1-n0) - g(x0,par).*n0))/(1+trapz(x0,n0)) )
fprintf('maximal force = %5.4f\n', trapz(x0,x0.*f(x0,par)./(f(x0,par)+g(x0,par)))  ) % Just to check that your final spring extension is correct
fprintf('force at tend = %5.4f\n', yout(numt,2*par.num+1))

end

%% ------------------------------------------------
% RHS of ode
function out = derivs(t,y,par)
x = y(1:par.num);
n = y(par.num+1:2*par.num);
L = y(2*par.num+1);

v = trapz(x, x.*(f(x,par).*(1-n) - g(x,par).*n))/(1+trapz(x,n));

out(1:par.num) = -v;
out(par.num+1:2*par.num) = (1-n).*f(x,par) - n.*g(x,par);
out(2*par.num+1) = v;
out = out';
end

%% ------------
function out = f(x,par)
out = 0 + (x>0 & x<=1).*par.f1.*x;
end

function out = g(x,par)
out = (x<0).*par.g2 + (x>0).*par.g1.*x;
end







