 % cirle map


 clear all
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)

% set parameters
gam = 0.6;
del = 1;
 delt=0.001;
t=[0:delt:6];
F = sin(pi*t).^4.*exp(gam*t);

G=F+del*exp(gam*t);

%do an iteration
 
tk(1) = 0.;
N=0;
j=1;
while N==0
    tj = tk(j);
Gkj=sin(pi*tj).^4.*exp(gam*tj)+del*exp(gam*tj);
Gk(j) = Gkj;
 if(isempty(find(F>=Gkj)))
     N=1;
 else
     j=j+1;
 ndx = min(find(F>=Gkj));
 tk(j) =t(ndx);
 end
end 
N=j;
gk(1) = 0;
for j = 1:N
    gk(2*j) = Gk(j);
    gk(2*j+1) = Gk(j);
    Tk(2*j-1) = tk(j);
    Tk(2*j) = tk(j);
end
figure(1)
plot(t,G,t,F,Tk,gk(1:2*N),'k--')
axis([0 4 0 20])
xlabel('t','fontsize',20)
text(3.65,8,'F(t)','fontsize',18)
text(3.65,16.5,'G(t)','fontsize',18)
text(0.3,3.75,'t_1','fontsize',18)

text(3.1,12.3,'t_3','fontsize',18)

text(1.3,5.95,'t_2','fontsize',18)
 
% plot some circle maps
gamlist=[0.75,0.694,0.67,0.55];
for ngam = 1:length(gamlist)
    clear t1 t2 t3 tk Gkj Tk
    gam=gamlist(ngam);

t=[0:delt:1];
F = sin(pi*t).^4.*exp(gam*t);
[M,I]=max(F);
ti=[0:delt:t(I)];
Gkj=sin(pi*ti).^4.*exp(gam*ti)+del*exp(gam*ti);
t1=[0:delt:3];
F1=sin(pi*t1).^4.*exp(gam*t1);
for j=1:I
    j1=min(find(F1>=Gkj(j)));  
     t2(j)= t1(j1);
end
t3=mod(t2,1);


% track a trajectory
tk(1)=0.35;
t=[0:delt:15];
F = sin(pi*t).^4.*exp(gam*t);
N=0;
j=1;
while N==0

    tj = tk(j);
Gkj=sin(pi*tj).^4.*exp(gam*tj)+del*exp(gam*tj);
 if(isempty(find(F>=Gkj)))
     N=1;
 else
     j=j+1;
 ndx = min(find(F>=Gkj));
 tk(j) =t(ndx);
 end
end
N=j;
 
for j = 1:N
    
    Tk(2*j-1) = tk(j);
    Tk(2*j) = tk(j);
end
Istop=max(find(t3(2:end)<t3(1:end-1)))
if(isempty(Istop))
figure(ngam+1)
plot(ti,t3,'r',t3,t3,'b--',mod(Tk(1:end-1),1),mod(Tk(2:end),1),'k--')
 axis([min(t3) max(t3) min(t3) max(t3)])
else
     
 figure(ngam+1)
plot(ti(1:Istop),t3(1:Istop),'r',ti(Istop+1:end),t3(Istop+1:end),'r',t3,t3,'b--',mod(Tk(1:end-1),1),mod(Tk(2:end),1),'k--')
 axis([min(t3) max(t3) min(t3) max(t3)]) 
end
xlabel('\Psi_n')
ylabel('\Psi_{n+1}')
box off
formatSpecF = '%6.2f\n';
 
 title(strcat('\gammaT = ',sprintf(formatSpecF,gam)),'fontsize',18)
end


 