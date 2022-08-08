clear all; close all; clc;

p.K1 = 10/1000; p.K2 = 10/1000; p.E2T = 50/1000;
p.v2 = p.E2T*1;

stimweight = 35/1000; 

% what is the response to a constant input?
steady_response = 1 - GK(stimweight,p)

% What is the response to an oscillatory input of the same total weight?
n = 250; 
response = zeros(1,n);
logfreq = linspace(-3,1,n);
freq = 10.^logfreq;
period = 1./freq;

for i=1:n
    p.T = period(i);
    p.w = 1;
    p.h = stimweight*p.T/p.w;
    h_comp(i) = ((1-GK(p.h,p))*p.w)/p.T;

    p.h = 0.35;
    p.w = stimweight*p.T/p.h;
    w_comp(i) = ((1-GK(p.h,p))*p.w)/p.T;
end
plot(logfreq,h_comp,'r',logfreq,w_comp,'m',logfreq,steady_response*ones(1,n),'b','LineWidth',2)
xlabel('frequency (1/T)')
ylabel('average response')
legend('changing amplitude','changing width','steady input')
set(gca,'FontSize',12)

%writematrix([freq' h_comp'],'adaptation.dat','Delimiter',' ')

function out=GK(H,p)
alpha = 1-H/p.v2;
beta = H*(p.K2+1)/p.v2 - (1-p.K2);
gamma = -p.K1;
out = (-beta + (beta.^2 - 4*alpha*gamma).^0.5)/(2*alpha);
end


