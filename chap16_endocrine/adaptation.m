clear all; close all; clc;

p.k1 = 0.5; p.k2 = 0.1;
p.km1 = 0.25; p.km2 = 0.1;

stimweight = 0.1; 

% what is the response to a constant input?
init = [0 0];
p.T = 1;
p.w = 1;
p.h = stimweight;
[t,y] = ode15s(@(t,y)rhs(t,y,p),[0,150],init);
steady_response = y(end,2);

% What is the response to an oscillatory input of the same total weight?
n = 2; 
response = zeros(1,n);
logfreq = linspace(-3,0,n);
freq = 10.^logfreq;
period = 1./freq;

for i=1:n
    p.T = period(i);
    p.h = 1;
    p.w = stimweight*p.T/p.h;
    options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',p.w/200);
    tspan = linspace(0,150*p.T,1000); % long enough to sit on the periodic solution
    init = [0 0];
    [t,y] = ode15s(@(t,y)rhs(t,y,p),tspan,init,options);
    %plot(t,y(:,2))
    %hold on
    init = [y(end,1) y(end,2)];
    tspan = linspace(t(end),t(end) + 2*p.T,1000);  % integrate over two more periods
    [t,y] = ode15s(@(t,y)rhs(t,y,p),tspan,init,options);
    %plot(t,y(:,2))
    response(i) = mean(y(:,2));
end
plot(logfreq,response,logfreq,steady_response*ones(1,n))

writematrix([freq' response'],'adaptation.dat','Delimiter',' ')


function out = rhs(t,y,p)
B = y(1);
D = y(2);
H = getH(t,p);
out(1) = p.k1*H*(1-B-2*D) - p.km1*B + 2*p.km2*D - 2*p.k2*(B^2);
out(2) = p.k2*(B^2) - p.km2*D;
out=out';
end

function out = getH(t,p)
out = 0;
if mod(t,p.T) < p.w
    out = p.h;
end
end