
%  -------------------------------------------------------------------
%
%   Reponse of the RyR model of Stern to a step increase in calcium.
%
%   For Chapter 7, Section 7.4.4 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

function RyR_Stern
%%   plotting of the steady state

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 18,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

global k1 km1 k2 km2
k1 = 35;
km1 = 0.06;
k2 = 0.5;
km2 = 0.01;

clist = [0.01:.01:0.5];
for j=1:length(clist)
    c=clist(j);

    [Rst(j),Ost(j),Ist(j)]=ROIss(c);

end

 
figure(2)
plot(clist,Rst,'r',clist,Ost,'b',clist,Ist,'g')
 hold on

% for external plotting
% keep = fplot(Ost,[0.01 0.5]);
% writematrix([keep.XData' keep.YData'],'test1.dat')

%% numerical solutions

%clear all

init = [1,0,0];                 % Assume that c = 0 initially. Probably not very accurate.
tspan = linspace(0,1.5,100);

figure(1)
clist = [1,2,3];
for j = 1:3
    c = clist(j);
 [T,Y] = ode15s(@(t,x)rhs(t,x,c),tspan,init);
plot(T,Y(:,2))
%writematrix( [T Y(:,2)],'test3.dat')   % for external plotting
hold on
end
box off
legend('boxoff')
legend('c=1 \mu M','c=2 \mu M','c=3 \mu M')
xlabel('time (ms)')
ylabel('O')
 hold off

% calculate the steady-state and peak responses

n = 50;
cc = linspace(0.0,0.5,n);
tspan = linspace(0,50,10000);
for i = 1:n
    c = cc(i);
    [T,Y] = ode15s(@(t,x)rhs(t,x,c),tspan,init);
    peak(i) = max(Y(:,2));
end

figure(2)
plot(cc,peak,'k')
xlabel('c (\mu M)')
ylabel('O')
  text(0.2,0.14,'steady response','fontsize',18)
text(0.25,0.75,'peak response','fontsize',18)

% for external plotting
% writematrix([cc',peak'],'test2.dat')

end % of main

%%
function out = rhs(t,x,c)
global k1 km1 k2 km2
R = x(1);
O = x(2);
I = x(3);
RI = 1 - R - O - I;

 out(1) = km1*O + km2*RI - R*(k1*c*c + k2*c);
out(2) = km2*I + k1*c*c*R - O*(k2*c + km1);
out(3) = k2*c*O + k1*c*c*RI - I*(km2 + km1);
out = out';

end


function[Rst,Ost,Ist]= ROIss(c)
global k1 km1 k2 km2
 
 % matrix equation
 A  =[[- km2-(k1*c*c + k2*c),km1- km2,- km2];[k1*c*c ,-(k2*c + km1),km2 ] ;[- k1*c*c,- k1*c*c+k2*c ,-(km2 + km1)- k1*c*c]]; 



Rhs = -[km2;0;  k1*c*c];

R =A\Rhs;
Rst=R(1);
Ost=R(2);
Ist=R(3);
end
