
% Hai-Murphy-Huxley model of smooth muscle. Use the method of
% characteristics to compute the distributions as functions of space and
% time.

function muscle_smooth
close all
clear all
clc

par.f1=200; 
par.g1=10;
par.g2=50;
par.h=1;

par.k1 = 2;
par.k2 = 1;
par.k5 = 2;
par.k6 = 1;
par.k7 = 0.1;

par.num = 250;      % number of space points
par.numt = 150;     % number of time outputs
par.tend = 5;       % final time

get_isometric_force(par);
get_time_dependent_force(par);

end

%%
function get_time_dependent_force(par)
% Use the method of characteristics to compute Nm, etc, as functions of x
% and t. Then calculate the force at each time. Starting with a relaxed
% muscle.

% initial conditions
x0 = linspace(-10,10,par.num);
point_spacing = 20/par.num;
par.xshift = floor(par.h/point_spacing) + 2;        % calculate how many points fit in length h and then add two more
par.delx = par.xshift*point_spacing;                % This is the distance between binding sites. Has to be bigger than h.

nm0 = ones(1,par.num);
nmp0 = zeros(1,par.num);
nam0 = zeros(1,par.num);
namp0 = zeros(1,par.num);
y0 = [x0 nm0 nmp0 nam0 namp0];

output_times = linspace(0,par.tend,par.numt);  % to make the curves for a constant v
% integrate the odes, one for each x point
[tout,yout]=ode15s(@(t,y)derivs_chars(t,y,par),output_times,y0);
for i=1:par.numt
    x = yout(i,1:par.num);
    nm = yout(i,par.num+1:2*par.num);
    nmp = yout(i,2*par.num+1:3*par.num);
    nam = yout(i,3*par.num+1:4*par.num);
    namp = yout(i,4*par.num+1:5*par.num);
    force(i) = trapz(x,x.*(nam+namp));
end
plot(tout,force)
hold on

% Now solve using the more complicated ODEs, with the sum on the RHS
[tout,yout]=ode15s(@(t,y)derivs_chars_alt(t,y,par),output_times,y0);
for i=1:par.numt
    x = yout(i,1:par.num);
    nm = yout(i,par.num+1:2*par.num);
    nmp = yout(i,2*par.num+1:3*par.num);
    nam = yout(i,3*par.num+1:4*par.num);
    namp = yout(i,4*par.num+1:5*par.num);
    force(i) = trapz(x,x.*(nam+namp));
end
plot(tout,force)

end

%% RHS of odes for the method of characteristics
function out = derivs_chars(t,y,par)
x = y(1:par.num);
nm = y(par.num+1:2*par.num);
nmp = y(2*par.num+1:3*par.num);
nam = y(3*par.num+1:4*par.num);
namp = y(4*par.num+1:5*par.num);

out(1:par.num) = -v(t);
out(par.num+1:2*par.num) = par.k2*nmp - par.k1*nm + par.k7*nam;                         % Nm
out(2*par.num+1:3*par.num) = par.k1*nm - (par.k2 + f(x,par)).*nmp + g(x,par).*namp;     % Nmp
out(3*par.num+1:4*par.num) = par.k6*namp - (par.k5+par.k7)*nam;                         % Nam
out(4*par.num+1:5*par.num) = par.k5*nam + f(x,par).*nmp - (par.k6+g(x,par)).*namp;      % Namp
out = out';
end

%% RHS of odes for the method of characteristics, but using the correct conservation law
function out = derivs_chars_alt(t,y,par)
x = y(1:par.num);
nm = y(par.num+1:2*par.num);
nmp = y(2*par.num+1:3*par.num);
nam = y(3*par.num+1:4*par.num);
namp = y(4*par.num+1:5*par.num);

% Now do all the shifts of nam and namp. By hand is the easiest to follow,
% although it's wordy
s=par.xshift; % to save typing
nams1 = zeros(par.num,1);  nams1(s+1:end) = nam(1:end-s);  
nams2 = zeros(par.num,1);  nams2(2*s+1:end) = nam(1:end-2*s);
nams3 = zeros(par.num,1);  nams3(3*s+1:end) = nam(1:end-3*s);
namsm1 = zeros(par.num,1);  namsm1(1:end-s) = nam(s+1:end);  
namsm2 = zeros(par.num,1);  namsm2(1:end-2*s) = nam(2*s+1:end);
namsm3 = zeros(par.num,1);  namsm3(1:end-3*s) = nam(3*s+1:end);

namps1 = zeros(par.num,1);  namps1(s+1:end) = namp(1:end-s);  
namps2 = zeros(par.num,1);  namps2(2*s+1:end) = namp(1:end-2*s);
namps3 = zeros(par.num,1);  namps3(3*s+1:end) = namp(1:end-3*s);
nampsm1 = zeros(par.num,1);  nampsm1(1:end-s) = namp(s+1:end);  
nampsm2 = zeros(par.num,1);  nampsm2(1:end-2*s) = namp(2*s+1:end);
nampsm3 = zeros(par.num,1);  nampsm3(1:end-3*s) = namp(3*s+1:end);

out(1:par.num) = -v(t);
out(par.num+1:2*par.num) = par.k2*nmp - par.k1*nm + par.k7*(nam + nams1 + nams2 + nams3 + namsm1 + namsm2 + namsm3);   % Nm
out(2*par.num+1:3*par.num) = par.k1*nm - (par.k2 + f(x,par)).*nmp + g(x,par).*(namp + namps1 + namps2 + namps3 + nampsm1 + nampsm2 + nampsm3);     % Nmp
out(3*par.num+1:4*par.num) = par.k6*namp - (par.k5+par.k7)*nam;                         % Nam
out(4*par.num+1:5*par.num) = par.k5*nam + f(x,par).*nmp - (par.k6+g(x,par)).*namp;      % Namp
out = out';
end

%%
function get_isometric_force(par) % Compute the isometric steady-state force. Just for fun.
x0 = linspace(-10,10,par.num);
k1=par.k1;
k2=par.k2;
k5=par.k5;
k6=par.k6;
k7=par.k7;
steady = (f(x0,par)*k1*(k5 + k7))./(f(x0,par)*k1*k5 + f(x0,par)*k1*k6 + f(x0,par)*k1*k7 + f(x0,par)*k6*k7 + g(x0,par)*k1*k5 + g(x0,par)*k2*k5 + g(x0,par)*k1*k7 + g(x0,par)*k2*k7 + k1*k6*k7 + k2*k6*k7) ...
            + (f(x0,par)*k1*k6)./(f(x0,par)*k1*k5 + f(x0,par)*k1*k6 + f(x0,par)*k1*k7 + f(x0,par)*k6*k7 + g(x0,par)*k1*k5 + g(x0,par)*k2*k5 + g(x0,par)*k1*k7 + g(x0,par)*k2*k7 + k1*k6*k7 + k2*k6*k7);
fprintf('isometric force is %5.4f\n',trapz(x0,x0.*steady))
end
%% ------------
function out = f(x,par)
out = 0 + (x>0 & x<=1).*par.f1.*x.*(1-x);
end

function out = g(x,par)
%out = (x<0).*par.g2 + (x>0).*par.g1.*x;
out = par.g2;
end

function out = v(t)
out = 0.05; %sin(10*t);
end







