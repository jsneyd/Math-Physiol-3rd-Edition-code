
function muscle_compute_max_v

close all
clear all
clc

par.f1=43.3; 
par.g1=10;
par.g2=209;
par.h=1;
par.v = 100; % constant velocity of contraction

%%  First calculate n for various v, and plot
xspan = linspace(2,-8,1500);
n0 = 0;
figure(1)
hold on
v = [1 10 50 100 120];
for i=1:5
par.v = v(i);
[xout,nout]=ode15s(@(x,n)derivs(x,n,par),xspan,n0);
plot(xout,nout,'LineWidth',2)
end
hold off
xlim([-2 2])
xlabel('x')
ylabel('n')
exportgraphics(gcf,'../../../Math-Physiol-3rd-edition/figures/chap_15_muscle/exercises/muscle_compute_max_v_1.pdf')

%% Next, plot the force-velocity curve for the given parameters
for i = 1:20
    v(i) = 0.1 + (i-1)*120/19;
    p(i) = get_load(v(i),par);
end
figure(2)
plot(p,v,'LineWidth',2)
xlabel('force')
ylabel('velocity')
exportgraphics(gcf,'../../../Math-Physiol-3rd-edition/figures/chap_15_muscle/exercises/muscle_compute_max_v_2.pdf')

%% Finally, plot the maximal velocity as a function of each parameter. Reset the parameter values for each plot.
par.f1=43.3; 
par.g1=10;
par.g2=209;
for i = 1:20
    par.g1 = 5 + (i-1)*10/19;
    g1(i) = par.g1;
    vmax(i) = fzero(@(v)get_load(v,par),100); % The maximal velocity occurs when the load = 0
end
figure(3)
plot(g1,vmax,'LineWidth',2)
xlabel('g_1')
ylabel('v_{max}')
exportgraphics(gcf,'../../../Math-Physiol-3rd-edition/figures/chap_15_muscle/exercises/muscle_compute_max_v_3.pdf')

par.f1=43.3; 
par.g1=10;
par.g2=209;
for i = 1:20
    par.g2 = 150 + (i-1)*100/19;
    g2(i) = par.g2;
    vmax(i) = fzero(@(v)get_load(v,par),100); % The maximal velocity occurs when the load = 0
end
figure(4)
plot(g2,vmax,'LineWidth',2)
xlabel('g_2')
ylabel('v_{max}')
exportgraphics(gcf,'../../../Math-Physiol-3rd-edition/figures/chap_15_muscle/exercises/muscle_compute_max_v_4.pdf')

par.f1=43.3; 
par.g1=10;
par.g2=209;
for i = 1:20
    par.f1 = 20 + (i-1)*70/19;
    f1(i) = par.f1;
    vmax(i) = fzero(@(v)get_load(v,par),100); % The maximal velocity occurs when the load = 0
end
figure(5)
plot(f1,vmax,'LineWidth',2)
xlabel('f_1')
ylabel('v_{max}')
exportgraphics(gcf,'../../../Math-Physiol-3rd-edition/figures/chap_15_muscle/exercises/muscle_compute_max_v_5.pdf')

end

%% -----------------------------------------------
function load = get_load(v,par)
% initial conditions
xspan = linspace(2,-8,1500);
n0 = 0;
par.v = v;
[xout,nout]=ode15s(@(x,n)derivs(x,n,par),xspan,n0);
load = -trapz(xout,xout.*nout);  % Has to be minus since x is going backwards in the integration
end

%% ------------------------------------------------
% RHS of ode
function out = derivs(x,n,par)
out = (-1/par.v)*( (1-n)*f(x,par) - n*g(x,par) );
end

%% ------------
function out = f(x,par)
out = 0 + (x>0 & x<1).*(par.f1*x);
end

function out = g(x,par)
out = par.g2*(x<=0) + par.g1*x.*(x>0);
end







