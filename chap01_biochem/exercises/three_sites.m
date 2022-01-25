clear all
close all
clc

syms c1 c2 c3 k1 km1pk2  k3 k4pkm3 k5 k6pkm5 e0  s  ee

fun1=k1*s*ee-(km1pk2)*c1-k3*s*c1+(k4pkm3)*c2==0;
fun2=k3*s*c1-(k4pkm3)*c2-k5*s*c2+(k6pkm5)*c3==0;
fun3=k5*s*c2-(k6pkm5)*c3==0;
ee=e0-c1-c2-c3;
                   
[c1,c2,c3]=solve([fun1,fun2,fun3],[c1,c2,c3])
