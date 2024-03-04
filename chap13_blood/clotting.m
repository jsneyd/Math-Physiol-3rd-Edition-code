
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
delt=0.1; delx=0.01; n=200; tend=350; nt=tend/delt;
K1=6.85; K2=11.0; K3=2.36; K4=0.087; K5=17.0; K6=0.066; D=2.6e-4;


% initial condition with stimulus on left
u1(1,1:n)=0; u1(1,1:10)=0.5;
u2(1,1:n)=0; u3(1,1:n)=0;

% use backward Euler integration

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
    plot(u1(1,:))
    hold on
    plot(u1(50/delt,:))
    plot(u1(250/delt,:),'red')
    plot(u1(300/delt,:),'green')
    plot(u1(350/delt,:),'black')
    box off
    hold off

%save('clotting')