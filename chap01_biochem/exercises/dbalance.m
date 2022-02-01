
clear all
close all
clc

syms a b k2 K1 k1 K2 C0 k3 K3  
fun1 = -(k2 + k1*K1)*a + k1*b + k2*K2*(C0 - a - b) == 0;
fun2 = -(k1 + k3*K3)*b + k1*K1*a + k3*(C0 - a - b) == 0;
[a,b] = solve([fun1,fun2], [a b]);

V=-a*k1*K1 + b*k1
[N,D]=numden(V)
