%the Goldwin negative feedback oscillator
function repressilator 

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
 'defaultlinelinewidth', 2.0);
global  alpha  p
%parameters
 
 
% The Golwin Oscillator
% For Keener and Sneyd, Mathematical Physiology, Third Edition, Chapter 10

% first get the HB curve
%stability for Goldwin model
P=[8.1:.01:20];

y=(8./(P-8)).^(1./P);

b = 1./((1+y.^P).*y);

 alpha = 0.2;
 
  p = 9;
figure(1)
plot(b,P,alpha,p,'*')
 


  m1=0.7 ;
  m2=m1;
  m3=0.8;

  init=[m1,m2,m3];
total=100;
tstep = 0.1;
tic
%specify the output points
tspan = [0:tstep:total];
[T,S] = ode23s(@deRHS,tspan, init, odeset('maxstep',10));  
toc

figure(2)
plot(T,S(:,1),T,S(:,2),T,S(:,3))
legend('boxoff')
legend('x_1','x_2','x_3')
box off
function s_prime=deRHS(t,sol) 
global  alpha  p 
x1=sol(1); 
x2=sol(2);
x3=sol(3);

x1p =  1/(1+x3^p)/alpha - x1;
x2p =   (x1-x2);
x3p =   (x2-x3); 


s_prime =[x1p;x2p;x3p];
 