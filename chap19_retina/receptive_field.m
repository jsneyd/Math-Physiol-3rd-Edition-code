% Code to solve exercises on receptive fields in Chapter 19

clear all
close all
clc

s1 = 3;
s2 = 1;
g1 = 3;
g2 = 1;
c = 1;

t = linspace(-3,3,100);
f = @(x) g1*s1*exp(-s1*s1.*x.*x)/pi^0.5 - g2*s2*exp(-s2*s2.*x.*x)/pi^0.5;   % The receptive field as a sum of Gaussians
transient = @(x,t) f(x).*(t-x/c).*exp(x/c-t);                               % transient response assuming a non-instantaneous response
R = (g1/2)*(1 + erf(s1*c*t)) - (g2/2)*(1 + erf(s2*c*t));                    % response to a semi-infinite bar moving at speed c

% Calculate moving bar response numerically. You can probably do this more
% efficiently. Oh dear.
for i=1:100
steady_f(i) = integral(f,-Inf,t(i));
transient_f(i) = integral(@(x)transient(x,t(i)),-Inf,t(i));
end

figure(1)
plot(t,R,'LineWidth',2)
hold on
%plot(t,steady_f,'r--','LineWidth',2)
plot(t,steady_f + transient_f,'r--','LineWidth',2)


%% Now calculate the response to a moving bar, width w

t = linspace(-3,8,300);
w = 5;
c=1;
R = (g1/2)*erf(s1*c*t) - (g2/2)*erf(s2*c*t) - (g1/2)*erf(s1*(c*t-w)) + (g2/2)*erf(s2*(c*t-w));                   

figure(2)
plot(t,R,'LineWidth',2)
set(gca,'Fontsize',14)
xlabel('t')
ylabel('r(t)')
box off
hold on

saveas(2,'../../Math-Physiol-3rd-Edition/figures/chap_19_retina/exercises/receptive_field.png')