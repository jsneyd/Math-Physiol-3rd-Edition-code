%  -------------------------------------------------------------------
%
%   Code to calculate and plot an effective diffusion coefficient.
%
%   For Chapter 8, Figure 8.21 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
% 
%   Written by James Keener and James Sneyd
% 
%  ------------------------------------------------------------------- 

clear all
close all
clc

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

del=0.1;
dellist = [0.2,0.1,0.01];
llist=[0.01:0.01:4];

for delj = 1:3
     del=dellist(delj);

    for j = 1:length(llist)
        l=llist(j);
       Deff(j) =  deff(l,del);
    end
    
    figure(1)
    plot(llist,Deff)
    hold on
    xlabel('l/L')
    ylabel('D_{eff}/D')
end
 
text(0.8,0.5,'\Delta=0.2','fontsize',18)
text(1.2,0.3,'\Delta=0.1','fontsize',18)  
text(2,0.135,'\Delta=0.01','fontsize',18)
hold off

clear Deff
% a different loop
dellist = [0:0.01:1];
llist = [0.000001,0.2,0.5];
for lj=1:3
    l = llist(lj);
    for dj = 1:length(dellist)
        del=dellist(dj);
        Deff(dj) =  deff(l,del);
    end
    if(lj==1)
        p=dellist.*(Deff-1);
        q=polyfit(dellist(5:end-1),p(5:end-1),1)
        pfit=polyval(q,dellist);
    
        figure(3)
        plot(dellist,p, dellist,pfit,'--')
    end

   figure(2)
   plot(dellist,Deff)
   hold on
   xlabel('\Delta')
   ylabel('D_{eff}/D')
  
end

text(0.1,0.97,'l/L=0','fontsize',18)
text(0.08,0.8,'l/L=0.2','fontsize',18)
text(0.18,0.5,'l/L=0.5','fontsize',18)
hold off
       
%%
function out=deff(l,del)
N=450;
for kj = 1:N
    for n = 1:N 
        k = 1+2*(kj-1);
        A(kj,n) = (2*n*tanh(2*n*pi*l*del)/(k*tanh(k*pi*l*(1-del)))+1)*n/(4*n^2-k^2);
    end
    R(kj) = 1/(2*pi*k^2);
    P(kj) = tanh(2*kj*pi*l*del);
end 
C=A\R';
out =P*C/l+del;
end
