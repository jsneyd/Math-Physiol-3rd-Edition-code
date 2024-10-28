
%    -------------------------------------------------------------------
%
%     Adaptation in a hormone receptor, shown by pulsatile stimulation 
%     at different frequencies.
%
%     For Chapter 16, Section 16.8.1 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------

function adaptation

clear all; close all; clc;

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
p.k1 = 0.5; p.k2 = 0.1;
p.km1 = 0.25; p.km2 = 0.1;

stimweight = 0.1;

% what is the response to a constant input?
% use quadratic formula to find steady state

p.T = 1;
p.w = 1;
p.h = stimweight;

alph = 2*p.k1*p.k2*p.h/p.km2;
bet = p.k1*p.h+p.km1;
gam = - p.k1*p.h;
disc = sqrt(bet^2-4*alph*gam);
root1 = (-bet+disc)/(2*alph);
steady_response = root1^2*p.k2/p.km2

% What is the response to an oscillatory input of the same total weight?
n = 40;
response = zeros(1,n);
logfreq = linspace(-3,0,n);
freq = 10.^logfreq;
period = 1./freq;
init = [0 0];

for i=1:n
    p.T = period(i);
    p.h = 1;
    p.w = stimweight*p.T/p.h;
    options = odeset('RelTol',1e-6,'AbsTol',1e-7,'MaxStep',p.w/200);
    tspan = linspace(0,30*p.T,1000); % long enough to sit on the periodic solution
    [t,y] = ode15s(@(t,y)rhs(t,y,p),tspan,init,options);

    %hold on
    init = [y(end,1),y(end,2)];
    tspan = linspace(t(end),t(end) + 2*p.T,1000);  % integrate over two more periods
    [t,y] = ode15s(@(t,y)rhsq(t,y,p),tspan,[init,0],options);
%     figure(1)
%         plot(y(:,1),y(:,2) ) %plot(t,y(:,2))
%         axis([0  0.6 0  0.2])
    response(i) = y(end,3)/(2*p.T);
end

figure(2)
    semilogx( freq,steady_response*ones(1,n),'--', freq,response)
    box off
    xlabel('Frequency of stimulation (1/T)')
    ylabel('Average response over period')
    legend('boxoff')
    legend('steady input','pulsatile input')
%writematrix([freq' response'],'adaptation.dat','Delimiter',' ')

end % of main

%%
function out = rhs(t,y,p)
    B  = y(1);
    D= y(2);

    H = getH(t,p);
    out(1) = p.k1*H*(1-B-2*D)-p.km1*B+2*p.km2*D-2*p.km2*B^2;
    out(2) = -2*p.km2*D+2*p.km2*B^2;
    out=out';
end

%%
function out = rhsq(t,y,p)
    % with quadrature
    B  = y(1);
    D = y(2);

    H = getH(t,p);
    out(1) = p.k1*H*(1-B-2*D)-p.km1*B+2*p.km2*D-2*p.km2*B^2;
    out(2) = -2*p.km2*D+2*p.km2*B^2;
    out(3) = D;

    out=out';
end

%%
function out = getH(t,p)
    out = 0;
    if mod(t,p.T) < p.w
        out = p.h;
    end
end
