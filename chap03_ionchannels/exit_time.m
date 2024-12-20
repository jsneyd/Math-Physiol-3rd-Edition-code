
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code to calculate k_1 and k_{-1} using mean first exit times
% for a double well potential

%   For Chapter 3,  Figures 3.18 A,B
%
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 
 

function exit_time
 
clear all
close all
clc
global A  

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

Ntot=100; 
dGlow=0.1; dGhigh=5;  % these are the limits for the amplitude A
deltaG = linspace(dGlow,dGhigh,Ntot);
b = 2; %location of minimum
 d=@(x)x;
 fun = @(x,y)exp(UU(x)).*exp(-UU(y));
 
for count=1:Ntot
    A=deltaG(count);  % go through a range of amplitudes
tau_1= quad2d(fun,0,1,-10,d);
k1(count) =1/tau_1;
     
    k1_approx(count) = ((-UUpp(0)*UUpp(1))^0.5/(pi))*exp(-A);


     
    tau_m1 = quad2d(fun,1,2,d,10);
    km1(count) = 1/tau_m1; 
 

    
    km1_approx(count) =  ((-UUpp (b )*UUpp(1))^0.5/(pi))*exp(-2*A);
end
 
 
figure(1)
      semilogy(deltaG,k1,'r',deltaG,k1_approx,'b--',deltaG,km1,'r',deltaG,km1_approx,'b--')
    xlabel('a=-\Delta G^0/k_BT')
    text(3,0.01,'k_{-1}','fontsize', 18)
    text(4.6,0.12,'k_1','fontsize', 18)
    box off
    fiddle=((-UUpp(b)*UUpp(1))^0.5)/((-UUpp(0)*UUpp(1))^0.5);
    
figure(3)
semilogy(deltaG,km1./k1,'r',deltaG,fiddle*exp(-deltaG),'b--')
   xlabel('a=-\Delta G^0/k_BT')
    ylabel('k_{-1}/k_1')
box off

dump = [deltaG' k1' k1_approx' km1' km1_approx' km1'./k1' fiddle*exp(-deltaG)'];
save('delG.dat','-ascii','dump')


Ntot=100; templow=0.03; temphigh=10;
A = 1;
temp=linspace(templow,temphigh,Ntot);
 
%dump = [temp' koff' koff_approx' kon' kon_approx'];
%save('delT.dat','-ascii','dump')
 
%%%%%%%%%%%%%%%%%%%%%%
function out=UU(x)
global A  
out = A*( 0.1319*(x.^6) - 0.0417*(x.^5) - 0.5347*(x.^4) - 1.3333*(x.^3) + 2.7778*(x.^2) );
%%%%%%%%%%%%%%%%%%%%%%
function out=UUpp(x)
global A  
out = A*( 6*5*0.1319*(x.^4) - 5*4*0.0417*(x.^3) - 4*3*0.5347*(x.^2) - 3*2*1.3333*(x.^1) + 2*2.7778 );
 