
%  -------------------------------------------------------------------
%
%   Program to solve the Llinas model of synaptic suppression.
%
%   For Chapter 8, Section 8.1.2 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd
%
%  -------------------------------------------------------------------

function Llinas

clear all
close all
clc

set(0,                       ...
'defaultaxesfontsize', 20,   ...
'defaultaxeslinewidth', 2.0, ...
'defaultlinelinewidth', 2.0);

global k10 k2 z1 s0 PCa F   ci ce n Vss k1  ohss tsw FbyRT

k10=2; k2=1; z1=1; s0=100; PCa=0.00001; F=96490; R=8.315; Tp=300;
ci=0.00001; ce=40; n=5;
Vss=-70;
FbyRT = F/(R*Tp*1000);  % units of mV^(-1)
k1=k10*exp(FbyRT*z1*Vss);
ohss= k1/(k1+k2);
tsw = 2.5;

t = [0:.01:5];
Vlist = [-40,-10,20,50,120];  % voltage in units of mV

figure(1)
for j = 1:length(Vlist)
    V = Vlist(j);
    oh = prob(t,V);
    currents= fCa(V,oh);
    ICa=currents(1,:);

    plot(t,-ICa,'LineWidth',2)
    xlabel('time (ms)')
    ylabel('-I_{Ca}(pA/(\mum)^2)')
    legend('boxoff')
    hold on
    ICa_save(j,:) = ICa;            % for external plotting
end
legend('V=-40mV','V=-10mV','V=20mV','V=50mV','V=120mV')

%save('Llinas.mat','t','ICa_save')   % for external plotting

x = linspace(-70,70,1000);  % units mV
phi_x=2*FbyRT*x;
Jss=PCa*2*F*phi_x.*(ci-ce*exp(-phi_x))./(1-exp(-phi_x));
Oss=k10*exp(FbyRT*z1*x).*(1./(k10*exp(FbyRT*z1*x)+k2));
ICass=Jss.*(s0/n).*((Oss).^n);

figure(2)
[hax, h1, h2] = plotyy(x,ICass,x,Oss);
ylabel(hax(1),'I_{Ca}')
ylim(hax(2),[0,1])
xlim([-70,70])
xlabel('V (mV)')
hold on
plot(x,Jss,'--');

% ylabel('\hat{o}')
legend('boxoff')
legend('I_{Ca}','j','\delta','Location','northwest','fontsize',18)


%save('Llinas.mat','x','ICass','Jss','Oss',"-append")  % for external plotting

for j = 1:length(Vlist)
    V = Vlist(j);
    Vt=(t<tsw).*V+(t>=tsw).*Vss;
    oh = probde(t,V);
    ICa = fCa(Vt,oh');

    ICa_save2(j,:) = ICa;       % for external plotting
    oh_save2(j,:) = oh;          % for external plotting

    figure(3)
    semilogy(t,-ICa,'LineWidth',2)
    xlabel('time (ms)')
    ylabel('-I_{Ca}(pA/(\mum)^2)')
    hold on
    legend('boxoff')

    figure(4)
    plot(t,oh)
    hold on
    legend('boxoff')

    xlabel('time (ms)')
    ylabel('open probability')
end

%save('Llinas.mat','ICa_save2','oh_save2',"-append")   % for external plotting


figure(3)
legend('V=-40mV','V=-10mV','V=20mV','V=50mV','V=120mV')
hold off
figure(4)
legend('V=-40mV','V=-10mV','V=20mV','V=50mV','V=120mV')
hold off


end % of main

%%     Calculate oh analytically
function oh=prob(t,V)
global k10 k2 z1  ohss FbyRT

k1=k10*exp(FbyRT*z1*V);
oh=k1./(k1+k2)*(1-exp(-(k1+k2)*t))+ohss*exp(-(k1+k2)*t);

end

%%
function currents=fCa(V,oh)
global  s0 PCa F ci ce n FbyRT

phi=2*FbyRT*V;
jj=PCa*2*F.*phi.*(ci-ce.*exp(-phi))./(1-exp(-phi));
ICa =(s0/n).*jj.*(oh.^n);

currents = [ICa];
end

%%  Calculate oh using the differential equation
function oh=probde(t,V)
global k10 k2 z1  Vss tsw FbyRT k2

k1=k10*exp(FbyRT*z1*Vss);
ohss= k1/(k1+k2);
sinit=ohss;
tspan = t;

[T,S] = ode15s(@(t,x)rhs(t,x,V),tspan,sinit);
% warning:  for piecewise continuous de's DO NOT use ode45 !!
oh =S(:,1);

end

%%
function out=rhs(t,s,V)
global Vss k10 k2 z1 tsw FbyRT

Vt=(t<tsw).*V+(t>=tsw).*Vss;
k1=k10*exp(FbyRT*z1*Vt);

out = k1.*(1-s)-k2*s;
end

