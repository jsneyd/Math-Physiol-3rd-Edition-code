
%    -------------------------------------------------------------------
%
% Plot the singular perturbation solution of the standing gradient osmotic
% solution.
%
%     For Chapter 18, Section 18.3.1 of
%     Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%     Written by James Keener and James Sneyd.
%
%    -------------------------------------------------------------------
 
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 2.0, ...
   'defaultpatchlinewidth', 0.7);


%parameters
D=1000;
r=0.05;
c0 = 0.3;
P=0.2;
L=100;
 
alp=0.1;

N0=0.3;

eps = D*r/(L^2*c0*P)
n0 = N0/(c0^2*P);
ud= 1/2 + sqrt(1 + 4*n0)/2;
wa = 2*alp*(ud-1);

% the outer solution
wo=wa*ud*[0:.001:1];
uo=wa*ud./wo;
y = (wa-wo-wa*ud*log((wo-wa*ud)./(wa-wa*ud)))/2+alp;
figure(1)
ndx = find(y>=alp&y<=1 );
w1=wo(ndx(end));
u1=uo(ndx(end));
tmp = (1-u1)*exp(w1*(y-1)/eps);
Y = [0,alp,y(ndx)]; 
Uo = [ud,ud,uo(ndx)+tmp(ndx)];
plot(y(ndx),uo(ndx) ,'--',Y,Uo)
legend('boxoff')
legend('outer solution','composite solution')
xlabel('y')
ylabel('u')
