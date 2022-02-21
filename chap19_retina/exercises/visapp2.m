%%
clear all
close all
clc

% This first bit calculates the inverse Fourier transform of the transfer
% function
syms a w trans

trans = 1/(a+i*w)^4;
latex(simplify(ifourier(trans)))
pretty(simplify(ifourier(trans)))


% This next bit plots the impulse responses. Just for convenience, all the
% parameters are set to 1.
t = linspace(0,15,200);
twostage = t.*exp(-t);
fourstage = (t.^3/factorial(3)).*exp(-t);
tenstage = (t.^9/factorial(9)).*exp(-t);

plot(t,twostage,t,fourstage,t,tenstage,'Linewidth',2)
legend('two stages','four stages','ten stages')
xlabel('t')
ylabel('K(t)')
set(gca,'FontSize',16)
box off