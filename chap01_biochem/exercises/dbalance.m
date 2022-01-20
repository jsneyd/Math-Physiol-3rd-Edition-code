
clear all
close all
clc

syms a b k2 km1 k1 km2 C0 k3 km3
fun1 = -(k2 + km1)*a + k1*b + km2*(C0 - a - b) == 0;
fun2 = -(k1 + km3)*b + km1*a + k3*(C0 - a - b) == 0;
[a,b] = solve([fun1,fun2], [a b]);

simplify(-a*km1 + b*k1)
