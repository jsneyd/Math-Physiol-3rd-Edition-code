% bistable waves  phase plane
%  -------------------------------------------------------------------
%
%  % find the traveling wave solution of the bistable equation using
%  shooting 
%
%   For Chapter 6, Figure 1.8, 1.9 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 

 clear all; close all; clc;
 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);  

 
global alf c
 
% 
cmykflg=1;
%if cmykflg=1, this code will make cmyk figure files
%cmykflg = 1;

alf = 0.1;
clist = [0;0.56;0.57;1;0.5665984]
for j = 1:5
c = clist(j);

lam = (-c+sqrt(c^2+4*alf))/2;
u0 = 0.001;
w0 = lam*u0;

s = [u0;w0];
tstep = 0.1;
t_end = 25;

tspan = [ 0:tstep:t_end];
[T,S] = ode23s(@deRHS,tspan,s,odeset('Events',@events)); 


figure(1)
if (j<5)
plot(S(:,1),S(:,2),'linewidth',2)
 hold on
 axis([0 1 0 0.2])
else
    plot(S(:,1),S(:,2),'--','linewidth',2)
 
 xlabel('U','fontsize',20)
 ylabel('W','fontsize',20)
end
text(0.08,0.18,'c=1','fontsize',20)
text(0.15,0.02,'c=0','fontsize',20)
text(0.67,0.02,'c=0.56','fontsize',20)
text(0.8,0.14,'c=0.57','fontsize',20)
 box off
end
hold off

 

[tmx,tj] = max(S(:,2));
tshift = T(tj);
figure(2)
plot(T-tshift,S(:,2),T-tshift,S(:,1),'linewidth',2)
text(3.5,0.1,'W(\zeta)','fontsize',20)
text(3.5,0.85,'U(\zeta)','fontsize',20)
 xlabel('\zeta','fontsize',20)
 axis([-10 10 0 1])
 box off
function s_prime=deRHS(t,s)
global alf c
v = s(1);
w = s(2);
 
Fv = w;
Fw = c*w-v*(1-v)*(v-alf);

s_prime = [Fv Fw]';
end

% event function
function [value,isterminal,direction]=events(t,s)
v = s(1);
w = s(2);

value = w*(1-v); 
isterminal = 1;
direction = 0;
end