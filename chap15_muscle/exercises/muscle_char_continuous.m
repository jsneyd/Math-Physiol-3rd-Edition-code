

function muscle_char
close all
clear all
clc

par.f1=43.3;
par.g1=10;
par.g2=209;
par.h=1;
par.v = 50;    % the velocity. Change the value to suit
par.num = 100; % number of space points


% initial conditions
x0 = linspace(-1,3,par.num);
n0 = f(x0,par)./(f(x0,par) + g(x0,par));
y0 = [x0 n0]; 

output_times = linspace(0,0.04,5);  % to make the curves for a constant v
% integrate the odes, one for each x point
[tout,yout]=ode45(@(t,y)derivs(t,y,par),output_times,y0);

% All the real work is done now. The rest is just plotting the output.
figure(1)
subplot(2,2,1)
plot(yout(1,1:par.num),yout(1,par.num+1:2*par.num),'LineWidth',2,'r')
xlabel('x')
ylabel('n')
title('t=0 s')
ylim([0 1])

subplot(2,2,2)
plot(yout(3,1:par.num),yout(3,par.num+1:2*par.num),'LineWidth',2,'b')
xlabel('x')
ylabel('n')
title('t=0.02 s')
ylim([0 1])

subplot(2,2,3)
plot(yout(4,1:par.num),yout(4,par.num+1:2*par.num),'LineWidth',2,'g')
xlabel('x')
ylabel('n')
title('t=0.03 s')
ylim([0 1])

subplot(2,2,4)
plot(yout(5,1:par.num),yout(5,par.num+1:2*par.num),'LineWidth',2,'m')
xlabel('x')
ylabel('n')
title('t=0.04 s')
ylim([0 1])


end


%% ------------------------------------------------
% RHS of ode
function out = derivs(t,y,par)
x = y(1:par.num);
n = y(par.num+1:2*par.num);

N = trapz(x,n);
out(1:par.num) = -par.v;
out(par.num+1:2*par.num) = (1-N).*f(x,par) - n.*g(x,par);

end

%% ------------
function out = f(x,par)
out = 0 + (x>0 & x<1).*(par.f1*x);
end

function out = g(x,par)
out = par.g2*(x<=0) + par.g1*x.*(x>0);
end







