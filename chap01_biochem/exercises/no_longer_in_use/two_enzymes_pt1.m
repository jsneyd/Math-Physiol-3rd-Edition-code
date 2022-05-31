

clear all
close all
clc

syms   l1  l2 th1 th2 s

 fun=(l1 + s )*(l2+s )/(s*(th1*s + th2));
 Fi=int(fun,s)
 K=subs(Fi,s,1)
  
