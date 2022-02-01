clear all
close all
clc

syms c1 c2 c3 k1 k2 k3 k4 k5 k6 e0 s ee K1 K2 K3

ee = e0-c1-c2-c3;
fun1=k1*s*ee-k1*K1*c1-k3*s*c1+k3*K2*c2==0;
fun2=k3*s*c1-k3*K2*c2-k5*s*c2+k5*K3*c3==0;
fun3=k5*s*c2-k5*K3*c3==0;
                
[c1,c2,c3]=solve([fun1,fun2,fun3],[c1,c2,c3])

V = k2*c1+k4*c2+k6*c3

[N,D] = numden(V)