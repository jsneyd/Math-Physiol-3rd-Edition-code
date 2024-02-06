% Fluid absorption as a function of lumenal Na concentration,
% in nondimensional variables.

%Reproducing   figure 18.3 in the gastrointestinal chapter.

clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

beta = 1;

gam = 10;
rho = 1;
u0 = 1;

ui = linspace(0,10,1000);
fi = (ui.^3)./(1+ui.^3) ;
ul = ui + beta*fi;

a = rho;
b = rho + ul;
c = ui - u0 + (1-gam)*beta*fi;
y = (-b + (b.^2 - 4*a*c).^0.5)./(2*a);

plot(ul,y,'LineWidth',2)
xlabel('lumenal Na^+ concentration, u_l')
ylabel('Flow rate, y')
set(gca,'FontSize',16)
box off