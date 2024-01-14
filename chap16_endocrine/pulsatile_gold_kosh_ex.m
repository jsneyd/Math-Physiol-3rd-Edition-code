% response of golbeter-Koshland to pulsatile input  
% h fixed, vary w
clear all; close all; clc;
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.2, ...
   'defaultlinelinewidth', 2.0); 

p.K1 = 1/100; p.K2 = 1/100; p.E2T = 5/100;
p.v2 = p.E2T*1;

stimweight = 35/1000; 

% what is the response to a constant input?
steady_response = 1 - GK(stimweight,p)

% What is the response to an oscillatory input of the same total weight?
n = 100; 
response = zeros(1,n);
logfreq = linspace(-2.4,0,n);
freq = 10.^logfreq;
period = 1./freq;

init = [0];

for i=1:n
    p.T = period(i);
    p.h = 0.35;  % hold w fixed
    p.w = stimweight*p.T/p.h;  % vary w
    options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',p.w/200);
    tspan = linspace(0,15*p.T,1000); % long enough to sit on the periodic solution
     
    [t,y] = ode15s(@(t,y)rhs(t,y,p),tspan,init,options);
 
    %plot(t,y(:,2))
    %hold on
    init = [y(end,1)];
    tspan = linspace(t(end),t(end) + 2*p.T,1000);  % integrate over two more periods
    [t,y] = ode15s(@(t,y)rhsq(t,y,p),tspan,[init,0],options);
    figure(1)
    plot(t,y(:,1))  
   % axis([0  0.6 0  0.2])
    response(i) = y(end,2)/(2*p.T);
end
 
figure(2)
semilogx( freq,steady_response*ones(1,n),'--', freq,response)
box off
xlabel('Frequency (1/T)')
ylabel('Average response ')
legend('boxoff')
legend('steady input','pulsatile input')
for i=1:n
    p.T = period(i);
    p.w = 1;
    p.h = stimweight*p.T/p.w;
    h_comp(i) = ((1-GK(p.h,p))*p.w)/p.T;

    p.h = 0.35;
    p.w = stimweight*p.T/p.h;
    w_comp(i) = ((1-GK(p.h,p))*p.w)/p.T;
end
%plot(logfreq,h_comp,'r',logfreq,w_comp,'m',logfreq,steady_response*ones(1,n),'b','LineWidth',2)
%xlabel('frequency (1/T)')
%ylabel('average response')
%legend('changing amplitude','changing width','steady input')
%set(gca,'FontSize',12)

%writematrix([freq' h_comp'],'adaptation.dat','Delimiter',' ')

function out=GK(H,p)
alpha = 1-H/p.v2;
beta = H*(p.K2+1)/p.v2 - (1-p.K2);
gamma = -p.K1;
out = (-beta + (beta.^2 - 4*alpha*gamma).^0.5)/(2*alpha);
end

function out = rhs(t,y,p)
s  = y(1);
sst = 1-s;
 
H = getH(t,p);
out(1) = p.v2*sst/(p.K2+sst) -H*s/(p.K1+s);
 
end

function out = rhsq(t,y,p)
% with quadrature
s  = y(1);
sst = 1-s;
 
H = getH(t,p);
out(1) = p.v2*sst/(p.K2+sst) -H*s/(p.K1+s);
 

 out(2) = sst;
 
out=out';
end
function out = getH(t,p)
out = 0;
if mod(t,p.T) < p.w
    out = p.h;
end
end
