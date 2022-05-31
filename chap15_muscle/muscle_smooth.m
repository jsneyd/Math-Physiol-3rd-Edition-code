
% Hai-Murphy-Huxley model of smooth muscle. Use the method of
% characteristics to compute the distributions as functions of space and
% time. Use two different versions of the model.

function muscle_smooth
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
par.numt = 500;     % number of time outputs
par.tend = 30;       % final time

% calculate the shifts needed in the complicated model version
point_spacing = 20/par.num;
par.xshift = floor(par.h/point_spacing) + 2;        % calculate how many points fit in length h and then add two more
par.delx = par.xshift*point_spacing;                % the distance of the shift

% initial conditions. use the conservation law to determine Nmp
x0 = linspace(-10,10,par.num);
nm0 = zeros(1,par.num);
nam0 = zeros(1,par.num);
namp0 = zeros(1,par.num);
y0 = [x0 nm0 nam0 namp0];
output_times = linspace(0,par.tend,par.numt);  % times to collect the output for plotting stuff

% First solve using the odes of the simplified model, one for each x point
[tout,yout]=ode15s(@(t,y)derivs_chars(t,y,par),output_times,y0);
for i=1:par.numt
    x = yout(i,1:par.num);
    nam = yout(i,2*par.num+1:3*par.num);
    namp = yout(i,3*par.num+1:4*par.num);
    force(i) = trapz(x,x.*(nam+namp));
end
figure(1)
plot(tout,force)
hold on

% Now solve using the more complicated ODEs, with the detailed conservation law
[tout,yout]=ode15s(@(t,y)derivs_chars_alt(t,y,par),output_times,y0);
for i=1:par.numt
    x = yout(i,1:par.num);
    nam = yout(i,2*par.num+1:3*par.num);
    namp = yout(i,3*par.num+1:4*par.num);
    force(i) = trapz(x,x.*(nam+namp));
end
plot(tout,force)
end

%% RHS of odes for the method of characteristics. Simple version.
function out = derivs_chars(t,y,par)
x = y(1:par.num);
nm = y(par.num+1:2*par.num);
nam = y(2*par.num+1:3*par.num);
namp = y(3*par.num+1:4*par.num);
nmp = ones(par.num,1) - nm - nam - namp;

out(1:par.num) = -v(t);
out(1*par.num+1:2*par.num) = par.k2*nmp - par.k1*nm + gL(x,par).*nam;                       % Nm
out(2*par.num+1:3*par.num) = par.k6*namp - (par.k5+gL(x,par)).*nam;                         % Nam
out(3*par.num+1:4*par.num) = par.k5*nam + f(x,par).*nmp - (par.k6+g(x,par)).*namp;      % Namp
out = out';
end

%% RHS of odes for the method of characteristics. Complex version using the correct conservation law
function out = derivs_chars_alt(t,y,par)
x = y(1:par.num);
nm = y(par.num+1:2*par.num);
nam = y(2*par.num+1:3*par.num);
namp = y(3*par.num+1:4*par.num);

% Now do all the shifts of Nam. By hand is the easiest to follow, although it's wordy
s=par.xshift; % to save typing
nams1 = zeros(par.num,1);  nams1(s+1:end) = namp(1:end-s);  
nams2 = zeros(par.num,1);  nams2(2*s+1:end) = namp(1:end-2*s);
nams3 = zeros(par.num,1);  nams3(3*s+1:end) = namp(1:end-3*s);
namsm1 = zeros(par.num,1);  namsm1(1:end-s) = namp(s+1:end);  
namsm2 = zeros(par.num,1);  namsm2(1:end-2*s) = namp(2*s+1:end);
namsm3 = zeros(par.num,1);  namsm3(1:end-3*s) = namp(3*s+1:end);

% Finally, determine Nmp by using the conservation law. 
nmp = ones(par.num,1) - nm - namp - nam - nams1 - nams2 - nams3 - namsm1 - namsm2 - namsm3;

out(1:par.num) = -v(t);
out(par.num+1:2*par.num) = par.k2*nmp - par.k1*nm + gL(x,par).*nam + ...
                                             gL(x+par.delx,par).*nams1 + gL(x+2*par.delx,par).*nams2 + gL(x+3*par.delx,par).*nams3 + ...
                                             gL(x-par.delx,par).*namsm1 + gL(x-2*par.delx,par).*namsm2 + gL(x-3*par.delx,par).*namsm3;   % Nm
out(2*par.num+1:3*par.num) = par.k6*namp - (par.k5+gL(x,par)).*nam;                         % Nam
out(3*par.num+1:4*par.num) = par.k5*nam + f(x,par).*nmp - (par.k6+g(x,par)).*namp;      % Namp
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

function out = v(t)
%out = 0.6; 
out = sin(2*t);
end







