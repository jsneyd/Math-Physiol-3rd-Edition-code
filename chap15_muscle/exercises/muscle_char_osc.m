
function muscle_char_osc

close all
clear all
clc

par.f1=43.3; 
par.g1=10;
par.g2=209;
par.h=1;
par.num = 200; % number of space points


% initial conditions
x0 = linspace(-1,4,par.num);
n0 = f(x0,par)./(f(x0,par) + g(x0,par));
y0 = [x0 n0];

nump = 250;
output_times = linspace(0,0.4,nump);  % to make the curves for a constant v
% integrate the odes, one for each x point
[tout,yout]=ode45(@(t,y)derivs(t,y,par),output_times,y0);

% All the real work is done now. The rest is just plotting the output.

figure(1)
subplot(2,2,1)
plot(yout(60,1:par.num),yout(60,par.num+1:2*par.num),'LineWidth',2,'r')
xlabel('x')
ylabel('n')
title('t = 0.1 s')
ylim([0 1])

subplot(2,2,2)
plot(yout(80,1:par.num),yout(80,par.num+1:2*par.num),'LineWidth',2,'b')
xlabel('x') 
ylabel('n')
title('t = 0.13 s')
ylim([0 1])

subplot(2,2,3)
plot(yout(100,1:par.num),yout(100,par.num+1:2*par.num),'LineWidth',2,'g')
xlabel('x')
ylabel('n')
title('t = 0.16 s')
ylim([0 1])

subplot(2,2,4)
plot(yout(120,1:par.num),yout(120,par.num+1:2*par.num),'LineWidth',2,'m')
xlabel('x')
ylabel('n')
title('t = 0.19 s')
ylim([0 1])


figure(2)
for i=1:nump
    tt(i,:) = tout(i)*ones(1,par.num); % has to be a better way to do this. I'm too lazy.
    plot3(yout(i,1:par.num),tt(i,:),yout(i,par.num+1:2*par.num))
    hold on
end
hold off
xlabel('x')
ylabel('t')
zlabel('n')

% In Matlab (but not in Octave) you can put all these curves into a surface by
% using Delaunay triangulation. However, the result is pretty awful. With a lot
% more effort you could make a nicer surface, but it's not really worth it. The
% pseudo-surface you get from plotting all the cross-sections is good enough for now.
% 
% % Put the data into long vectors for use by the delaunay triangulation
% xx = reshape(yout(:,1:par.num),1,[]);
% nn = reshape(yout(:,par.num+1:2*par.num),1,[]);
% tt = reshape(tt,1,[]);
% 
% tri = delaunay(xx,tt);
% figure(2)
% trisurf(tri, xx, tt, nn, 'LineStyle','none');

end

%% ------------------------------------------------
% RHS of ode
function out = derivs(t,y,par)
x = y(1:par.num);
n = y(par.num+1:2*par.num);

out(1:par.num) = -v(t);
out(par.num+1:2*par.num) = (1-n).*f(x,par) - n.*g(x,par);
out = out';
end

%% ------------
function out = f(x,par)
out = 0 + (x>0 & x<1).*(par.f1*x);
end

function out = g(x,par)
out = par.g2*(x<=0) + par.g1*x.*(x>0);
end

function out = v(t)
out = 50*sin(50*t);  
end





