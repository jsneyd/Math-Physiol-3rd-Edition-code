%%
clear all
close all
clc

syms k w t x real
assume(k>0);
assume(x>0);

n=2; k=1;
Phat = k^n/((k+i*w)^n)
P = ifourier(Phat)
pretty(P)

% Just for fun, make sure that the integral of P (from 0 to infinity) is 1.
% Which it is, of course.

totalint = int(P,0,inf)

% Now calculate the mean of P
m = int(x*P,0,inf)                  % This comes out to be n/k
v = int((x-m)^2*P,0,inf)            % This comes out to be n/k^2
