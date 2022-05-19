

function muscle_char
close all
clear all
clc

% for constant stretch change sign of v as required
par.v = 50;

% v0=25;
% om=50;
% par.v=@(t) v0*sin(om*t);


% binding functions
par.f1=43.3;
par.g1=10;
par.g2=209;
par.h=1;
par.num = 100; % number of space points


% initial conditions
x0 = linspace(-1,3,par.num);
n0 = f(x0,par)./(f(x0,par) + g(x0,par));
y0 = [x0 n0];

output_times = linspace(0,0.04,5);
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

out(1:par.num) = -par.v;
out(par.num+1:2*par.num) = (1-n).*f(x,par) - n.*g(x,par);

end

%% ------------
function out = f(x,par)
out = 0 + (x>0 & x<1).*(par.f1*x);
end

function out = g(x,par)
out = par.g2*(x<=0) + par.g1*x.*(x>0);
end






