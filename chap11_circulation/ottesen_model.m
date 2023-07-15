close all
clear all
clc

ca=1.55; R=1.05; Vs=67.9; alpha=0.84; beta=1.17; mu=93;

Pbar=mu*((alpha/beta)^(1/7));
gprime_s = -7*(Pbar^6)*mu^7/((mu^7+Pbar^7))^2;
a=sqrt(-alpha*gprime_s*Vs);

gprime_p=-gprime_s;
b=sqrt(beta*gprime_p*Vs);

coeff1=ca*ca; coeff2=1/(R*R) - 2*ca*b*b; coeff3=b^4-a^4;

xi_up=(-coeff2 + sqrt(coeff2^2-4*coeff1*coeff3))/(2*coeff1);
xi_down=(-coeff2 - sqrt(coeff2^2-4*coeff1*coeff3))/(2*coeff1);

z_up=sqrt(xi_up)
z_down=sqrt(xi_down)
upper_bound = a*a*R

z=z_up;
tau=(asin(z/(a*a*R)))/z
tau=(acos((ca*z*z-b*b)/(a*a)))/z
