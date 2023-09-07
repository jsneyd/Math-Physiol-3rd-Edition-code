% hindmarsh Rose 3 variable polynomial system subsystem

function hindmarsh_rose 
global   Iapp  r s  x1
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);

%parameters
 
Iapp=0;
r=0.001;
s=4;
x1 = -(1+sqrt(5))/2; ...
    
 Ilist=[0.4,2,4];
 for j=1:length(Ilist)
     Iapp=Ilist(j)
%initial data
zstart=[0.2,1.8,2];

init=[0,0,zstart(j)];

total=2000;
tstep = 0.01;
tic
%specify the output points
tspan = [0:tstep:total];
[T,S] = ode23(@deRHS,tspan, init, odeset('maxstep',10));  
toc

figure(2*j-1)
plot(T ,S(:,1),T,S(:,3) ,'--')
  formatSpecF = '%6.2f\n';
 
 title(strcat('I_{app} = ',sprintf(formatSpecF,Iapp)),'fontsize',18)
legend('boxoff')
legend('x','z')
ylabel('x')
xlabel('t')

figure(2*j)
plot(S(:,3),S(:,1) )
 formatSpecF = '%6.2f\n';
 
 title(strcat('I_{app} = ',sprintf(formatSpecF,Iapp)),'fontsize',18)
xlabel('z')
ylabel('x')
 end
function s_prime=deRHS(t,sol) 
global  Iapp  r s x1
x=sol(1); 
y=sol(2);
z=sol(3);

xp = y - x^3 + 3*x^2 + Iapp -z  ;
yp = 1 - 5*x^2 - y;
zp = r*(s*(x-x1)-z);

s_prime =[xp;yp;zp];
