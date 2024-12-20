%  -------------------------------------------------------------------
%
%   The fast subsystem of the Hindmarsh-Rose model of electrical bursting.
%
%   For Chapter 9, Section 9.1.3 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd
%
%  -------------------------------------------------------------------

function hindmarsh_rose_fast

close all
clear all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

global Iapp

%parameters
% change Iapp to see different phase portraits:
 Iapplist=[0,0.4,2 ];

 for j = 1:3
     Iapp=Iapplist(j);

 % make a phase plane
x  = [-2:.01:2];
y1 = x.^3 - 3*x.^2 - Iapp ;
y2 = 1 - 5*x.^2;
x1 = -1/2-sqrt(5)/2;
x2 = -1/2+sqrt(5)/2;
ys1 = 1 - 5*x1.^2;
ys2 = 1 - 5*x2.^2;

% initial data etc
init=[0.62,-0.86];
total=100;
tstep = 0.01;
tspan = [0:tstep:total];
[T,S] = ode23(@deRHS,tspan, init, odeset('maxstep',10));

 
figure(j)
plot(S(:,1),S(:,2),x,y1,'--',x,y2,'--')
axis([-2 2 -8 2])
box off 
legend('boxoff')
legend('solution trajectory','dx/dt=0','dy/dt=0','location','northwest')
ylabel('y')
xlabel('x')
formatSpecF = '%6.2f\n';
     
    title(strcat('I_{app} = ',sprintf(formatSpecF,Iapp)),'fontsize',18)
if (Iapp==0)
    hold on
    text(x1,ys1,'*','fontsize', 25)
    text(x2,ys2,'*','fontsize', 25)
    hold off
end
 end
 end % of main


%%
function s_prime=deRHS(t,s)
global Iapp
x=s(1);
y=s(2);

xp = y - x^3 + 3*x^2 + Iapp ;
yp = 1 - 5*x^2 - y;

s_prime =[xp;yp];
end
