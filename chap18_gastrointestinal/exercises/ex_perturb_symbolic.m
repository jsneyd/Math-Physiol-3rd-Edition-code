clear all
close all
clc

syms m0 a y f1 f2

f1 = -m0*y*y + 1 + m0*a*(2-a);
f2 = -2*m0*a*y + 1 + 2*m0*a;

simplify(int(f1,y,[0 a]) + int(f2,y,[a 1]))