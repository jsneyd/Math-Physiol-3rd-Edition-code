
%    -------------------------------------------------------------------
%  
%     Hai-Murphy-Huxley model of smooth muscle contracting against a
%     linear spring.
%  
%     For Chapter 15, Section 15.9.2 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%  
%     Written by James Keener and James Sneyd.
%  
%    -------------------------------------------------------------------

close all
clear all
clc

global num k2 k6 f1 g1 g2 gL1 gL2

k2 = 0.1;
k6 = 0.1;


f1 = 0.88; 
g1 = 0.21;
g2 = 4.4;
gL1 = 0.01;
gL2 = 0.2;

num = 200;       % number of space points
numt = 100;      % number of time outputs
tend = 10;       % final time


% Initial conditions
x0 = linspace(-2, 2, num);
nm0 = ones(1, num);
nam0 = zeros(1, num);
namp0 = zeros(1, num);
y0 = [x0, nm0, nam0, namp0];  % the final ode is for the spring length
tout = linspace(0, tend, numt);

% solve the odes
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, soln] = ode15s(@deriv, [0 tend], y0, options);

% calculate and plot the force
figure(1)
p = zeros(1, length(t));
for i = 1:length(t)
    x = soln(i, 1:num);
    nam = soln(i, 2*num+1:3*num);
    namp = soln(i, 3*num+1:4*num);
    p(i) = trapz(x, x .* (nam + namp));
end
plot(t, p)

% save the results
% Uncomment the following line if you want to save the data
% save('test.dat', 't', 'p', '-ascii')

%%
function dydt = deriv(t, y)
    global num k2 k6

    x = y(1:num);
    nm = y(num+1:2*num);
    nam = y(2*num+1:3*num);
    namp = y(3*num+1:4*num);
    nmp = ones(num,1) - nm - nam - namp;

    ff = arrayfun(@f, x);
    gg = arrayfun(@g, x);
    ggL = arrayfun(@gL, x);

    dum = trapz(x, x .* (ff .* nmp - gg .* namp - ggL .* nam));
    v = dum / (1 + trapz(x, namp + nam));
    dydt(1:num) = -v;
    dydt(num+1:2*num) = k2 * nmp - k1(t) * nm + ggL .* nam;
    dydt(2*num+1:3*num) = k6 * namp - (k5(t) + ggL) .* nam;
    dydt(3*num+1:4*num) = k5(t) * nam + ff .* nmp - (k6 + gg) .* namp;
    dydt = dydt(:);
end

function val = f(x)
    global f1
    val = 0;
    if x > 0 && x <= 1
        val = f1 * x;
    end
end

function val = g(x)
    global g1 g2
    if x <= 0
        val = g2;
    else
        val = g1 * x;
    end
end

function val = gL(x)
    global gL1 gL2
    if x <= 0
        val = gL2;
    else
        val = gL1 * x;
    end
end

function val = k1(t)
    val = 0.35 - 0.34 * sin(2 * t);
end

function val = k5(t)
    val = 0.35 - 0.34 * sin(2 * t);
end


