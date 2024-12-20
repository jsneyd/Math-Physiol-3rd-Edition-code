
% ------------------------------------------------
% Responses of receptive fields to moving bars and steps.
%
% For Chapter 19, Section 19.6, as well as a number of exercises.
%
% Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% Written by James Keener and James Sneyd
%
% ------------------------------------------------

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);

s1 = 3;
s2 = 1;
g1 = 3;
g2 = 1;
c = 1;
global w
t = linspace(-3,3,100);
f = @(x) g1*s1*exp(-s1*s1.*x.*x)/pi^0.5 - g2*s2*exp(-s2*s2.*x.*x)/pi^0.5;   % The receptive field as a sum of Gaussians
transient = @(x,t) f(x).*(t-x/c).*exp(x/c-t);                               % transient response assuming a non-instantaneous response
R = (g1/2)*(1 + erf(s1*c*t)) - (g2/2)*(1 + erf(s2*c*t));                    % response to a semi-infinite bar moving at speed c

% Calculate moving bar response numerically. You can probably do this more
% efficiently.  
for i=1:100
    steady_f(i) = integral(f,-Inf,c*t(i));
    transient_f(i) = integral(@(x)transient(x,t(i)),-Inf,c*t(i));
end

figure(1)
    plot(t,R,'LineWidth',2)
    hold on
    %plot(t,steady_f,'r--','LineWidth',2)
    plot(t,steady_f + transient_f,'r--','LineWidth',2)
    hold off

%% Now calculate the response to a moving bar, width w

t = linspace(-3,8,300);
w = 5;
clist = [1,10];
for j = 1:2
    c=clist(j);
    R =  (g1/2)*erf(s1*c*t) - (g2/2)*erf(s2*c*t) - (g1/2)*erf(s1*(c*t-w)) + (g2/2)*erf(s2*(c*t-w)); 
    figure(2)
        plot(t,R, 'LineWidth',2)
        set(gca,'Fontsize',14)
        xlabel('t')
        ylabel('r(t)')
        box off
        hold on
end

legend('c=1','c=10')
hold off

%saveas(2,'../../Math-Physiol-3rd-Edition/figures/chap_19_retina/exercises/receptive_field.png')

%% Now calculate the response to a step and a bar, when the response is not instantaneous

c = 1;
w = 8;
n = 200;
t = linspace(-2,20,n);
imp = @(t)(1.2*t.*exp(-t) - 0.45*t.^2.*exp(-t)).*heaviside(t);

for i = 1:n
    step_response(i) = integral(@(s)imp(t(i)-s),0,t(i));
end

figure(3)
    plot(t,imp(t),t,step_response,'LineWidth',2)
    set(gca,'Fontsize',14)
    xlabel('t')
    ylabel('r(t)')
    box off
    legend('impulse response','step response')
%saveas(3,'../../Math-Physiol-3rd-Edition/figures/chap_19_retina/exercises/receptive_field_2.png')


%% Finally, calculate the response to a moving step, when the response is not instantaneous

% You have to be careful with the resolution (the parameter n) here. It's a tradeoff between speed and accuracy. 
% Octave does the integral a lot slower than Matlab does, and struggles with n=200, which looks a lot better than n=50.

w = 8;
n = 50;    
t = linspace(-2,20,n);
clist=[1,0.5,3];
for j = 1:3
    c=clist(j); 
    stim = @(x,t)(heaviside(t-x/c) - heaviside(t-x/c-w/c));
    for i = 1:n
        bar_response(i) = integral2(@(s,x)stim(x,t(i)).*imp(t(i)-s).*f(x),-Inf,t(i),-Inf,Inf);
    end   
    figure(4)
        plot(t,bar_response,'LineWidth',2)
        set(gca,'Fontsize',18)
        xlabel('t')
        ylabel('r(t)')
        box off
        hold on
end
legend('c=1','c=0.5','c=3')
hold off
%saveas(4,'../../Math-Physiol-3rd-Edition/figures/chap_19_retina/exercises/receptive_field_3.png')

 