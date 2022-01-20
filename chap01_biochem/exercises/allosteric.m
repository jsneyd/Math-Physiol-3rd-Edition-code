clear all
close all
clc

syms k1 km1 k2 km2 k3 km3 e0 i s x y z ss dd nn

fun1 = (k1*s + k3*i)*(e0-x-y-z) - km1*x - km3*y-k2*x==0;
fun2 = (km3+k1*s)*y - km1*z - k3*i*(e0-x-y-z) ==0;
fun3 = (km3 + km1)*z - k3*i*x - k1*s*y == 0;

[x,y,z] = solve([fun1,fun2,fun3],[x,y,z]);

[nn,dd] = numden(x);
ss = nn/factor(dd);
simplify(ss)