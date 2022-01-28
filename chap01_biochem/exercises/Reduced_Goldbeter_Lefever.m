clear all
close all
clc

syms c0 c1 c2 t0 r0 r1 r2 s1 s2  k k1 km1 k2 km2  k3 km3 k4 km4 E0
 
eqt =k1*r0-km1*t0==0;  
 
eqr0= -k1*r0+km1*t0 -k2*s1*r0+km2*c0 +k*c0 -2*k3*s2*r0+km3*r2 
eqr1= -k2*s1*r1+km2*c1 +k*c1+2*k3*s2*r0-km3*r1 -k3*s2*r1+2*km3*r2 
eqr2= -k2*s1*r2+km2*c2 +k*c2 +k3*s2*r1-2*km3*r2 
eqc0=k2*s1*r0-km2*c0 -k*c0 
eqc1=k2*s1*r1-km2*c1 -k*c1 
eqc2=k2*s1*r2-km2*c2 -k*c2 
        
%check conservation
                      
simplify(eqt+eqr0+eqr1+eqr2+eqc0+eqc1+eqc2);
                               0
eqT=t0+r0+r1+r2+c0+c1+c2-E0 ;
 
                   
[c0,c1,c2,t0,r0,r1,r2]=solve([eqt ,eqr0 ,eqr1 ,eqr2 ,eqc0 ,eqc1 ,eqT ],[c0,c1,c2,t0,r0,r1,r2])

eqs1=-k2*s1*r0+km2*c0-k2*s1*r1+km2*c1-k2*s1*r2+km2*c2
eqs2=-2*k3*s2*r0 -k3*s2*r1+km3*r1+2*km3*r2 +k*c0+k*c1 +k*c2 
simplify(eqs1+eqs2)
