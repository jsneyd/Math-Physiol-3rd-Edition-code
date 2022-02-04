% Fluid absorption as a function of lumenal Na concentration,
% in nondimensional variables.


clear all
close all
clc

% First reproduce the figure in the gastrointestinal chapter.

beta = 1;
gam = 10;
rho = 1;
u0 = 1;

ui = linspace(0,8,1000);
fi = (ui.^3)./(1+ui.^3) ;
ul = ui + beta*fi;

a = rho;
b = rho + ul;
c = ui - u0 + (1-gam)*beta*fi;
y = (-b + (b.^2 - 4*a*c).^0.5)./(2*a);

plot(ul,y,'LineWidth',2)
xlabel('lumenal sodium concentration, u_l')
ylabel('flow rate, y')
set(gca,'FontSize',14,'LineWidth',1.5)
box off
hold on

% Next add to the figure the curve assuming that the Na flux depends on the
## voltage

V = -1;
ul = ui*exp(V).*exp(beta*fi);

a = rho;
b = rho + ul;
c = ui - u0 + (1-gam)*beta*fi;
y = (-b + (b.^2 - 4*a*c).^0.5)./(2*a);

plot(ul,y,'LineWidth',2)

% Save the file. For convenience. You probably don't want this line.
saveas(1,'../../../Math-Physiol-3rd-Edition/figures/chap_18_gastrointestinal/exercises/ex_sodium_flux.png')
