% Code to simulate the direction detection model
clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);

%parameters
Dx = 0.5;
Dt = 0.2;
c = 0.2;

t=[-20:0.1:20];
cval = [0.2,-0.2];
for j = 1:2
    c = cval(j);
R = S(-c*t).*S(Dx-c*t+c*Dt)-S(Dx-c*t).*S(-c*t+c*Dt);
figure(1) 
plot(t,R)
hold on
xlabel('t')
ylabel('R(t)')
end
legend('boxoff')
legend('c=0.2','c=-0.2','location','northwest')
box off

hold off


function out =S(x)

out = 10*(1+exp(-(x-1).^2) + 0.8*exp(-(x+0.2).^2));
end

