
%  -------------------------------------------------------------------
%
%   Compute a traveling wave in the clotting model.
%
%   For Chapter 13, Section 13.4.3 of
%   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
%
%   Written by James Keener and James Sneyd.
%
%  -------------------------------------------------------------------

clear all
close all
clc
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0)
delt=0.1; delx=0.01; n=120; tend=250; nt=tend/delt;
K1=6.85; K2=11.0; K3=2.36; K4=0.087; K5=17.0; K6=0.066; D=2.6e-4;


% initial condition with stimulus on left
u1(1,1:n)=0; u1(1,1:10)=0.5;
u2(1,1:n)=0; u3(1,1:n)=0;

% use backward Euler integration
x = [1:n]*delx;
% matrices for the backward Euler
lam=D*delt/(delx*delx);
A=(1+2*lam)*diag(ones(1,n)) - lam*diag(ones(1,n-1),1) - lam*diag(ones(1,n-1),-1);
A(1,2)=-2*lam; A(n,n-1)=-2*lam;

for i=2:nt
    uu1=u1(i-1,:); uu2=u2(i-1,:); uu3=u3(i-1,:);
    rhs1= K1*uu1.*uu2.*(1-uu1).*(1+K2*uu1)./(1+K3*uu3) - uu1;
    rhs2= uu1-K4*uu2;
    rhs3=K5*uu1.*uu1 - K6*uu3;
    
    u1(i,:) = A\(uu1' + delt*rhs1');
    u2(i,:) = A\(uu2' + delt*rhs2');
    u3(i,:) = A\(uu3' + delt*rhs3');
end

%surface(u1)
figure(1)
   % plot(u1(1,:))
   
    plot(x,u1(250 ,:),'red')
    hold on
    plot(x,u1(1000 ,:),'red')
    ylabel('u_1')
   % plot(u1(300/delt,:),'green')
    %plot(u1(350/delt,:),'black')
    box off
     text(0.35,0.15,'t=25','fontsize',18)
      text(0.75,0.2,'t=100','fontsize',18)
    yyaxis right
     %plot(u3(1,:))
    hold on
    plot(x,u3(250 ,:),'--blue')
    plot(x,u3(1000 ,:),'--blue')
    %plot(u3(300/delt,:),'green')
    %plot(u3(350/delt,:),'black')
    hold off
    xlabel('x')
    ylabel('u_3')
   

    figure(2)
     plot(x,u1(2000 ,:),'red')
     ylabel('u_1')
      text(0.2,0.37,'t=200','fontsize',18)
    hold on
     yyaxis right
     plot(x,u3(2000,:),'--blue')
     xlabel('x')
     ylabel('u_3')
     legend('boxoff')
     legend('u_1','u_3')
     box off
%save('clotting')