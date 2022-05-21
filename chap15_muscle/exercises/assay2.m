clear all
close all
clc

syms A g2
% f1 = 65;
% g1 = 15;
% g2 = 10;

out1 = g2/2 + A/24;
out2 = g2/sqrt(2) - g2*g2/(2*A);

dd = out1-out2;
root = solve(diff(dd,g2)==0,g2)
simplify(subs(dd,g2,root))