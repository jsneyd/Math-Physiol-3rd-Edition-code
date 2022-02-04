% Fluid absorption as a function of lumenal Na concentration,
% in nondimensional variables.

% Just reproducing a figure in the gastrointestinal chapter.

clear all
close all
clc

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
xlabel('lumenal sodium concentration, u_l')
ylabel('flow rate, y')
set(gca,'FontSize',14)
box off