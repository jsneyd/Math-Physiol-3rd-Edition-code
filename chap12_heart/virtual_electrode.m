% this is the calculation to 

%  -------------------------------------------------------------------
%
%   Find virtual electrodes in a bidomain model
%
%   For Chapter 12, Section 12.4.3 of
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
   'defaultaxeslinewidth', 1.5, ...
   'defaultlinelinewidth', 2.0)
N=101; % number of grid points (must be odd)
Nsq = N^2;
dx = 14/N;

%Here are a number of cases to consider 
%conductances; anisotropy ratios are unequal
six = 1;
siy = 0.09;
sex = 0.35;
sey = 0.11;


%equal anisotropy
six = .1;
sex = .2;
siy = 1;
sey = 2;

% inverse anisotropy
six = 2;
siy = .2;
sex = .2;
sey = 2;

% nominal values of Roth
six = 2;
siy = .2;
sex = 8;
sey = 2;

%my values
six = 1;
siy = .1;
sex = .35;
sey = 1.33;

%corrected values
six = 1;
siy = .05;
sex = .35;
sey = 1.33;


offdiag1 = [ones(N-2,1);2]; % off-diagonal -1
offdiag2 = [2;ones(N-2,1)]; % off-diagonal +1
for j = 1:N-1
    offdiag1 = [offdiag1;0;ones(N-2,1);2];
    offdiag2 = [offdiag2;0;2;ones(N-2,1)];
end
spoffdiag1 = [offdiag1;0];
spoffdiag2 = [0;offdiag2];

offdiag3 = [ones(N*(N-2),1);2*ones(N,1)]; 
spoffdiag3 = offdiag3;
offdiag4 = [2*ones(N,1);ones(N*(N-2),1)];
spoffdiag4 = [zeros(N,1);offdiag4];


A = (-2*(six+siy)*spdiags(ones(Nsq,1),0,Nsq,Nsq)+six*spdiags(spoffdiag2,1,Nsq,Nsq) ...
    + six*spdiags(spoffdiag1,-1,Nsq,Nsq) + siy*spdiags(spoffdiag3,-N,Nsq,Nsq) ...
+siy*spdiags(spoffdiag4,N,Nsq,Nsq))/dx^2;

B = (-2*(sex+sey)*spdiags(ones(Nsq,1),0,Nsq,Nsq)+sex*spdiags(spoffdiag2,1,Nsq,Nsq) ...
    + sex*spdiags(spoffdiag1,-1,Nsq,Nsq)+ sey*spdiags(spoffdiag3,-N,Nsq,Nsq) ...
    +sey*spdiags(spoffdiag4,N,Nsq,Nsq))/dx^2;

C = spdiags(ones(Nsq,1),0,Nsq,Nsq);

M = [A-C,C;C,B-C];

%AA = (-2*(six+siy)*diag(ones(Nsq,1))+six*diag(offdiag2,1) ...
%    + six*diag(offdiag1,-1) + siy*diag(offdiag3,-N) ...
%+siy*diag(offdiag4,N))/dx^2;

%BB = (-2*(sex+sey)*diag(ones(Nsq,1))+sex*diag(offdiag2,1) ...
%    + sex*diag(offdiag1,-1)+ sey*diag(offdiag3,-N) ...
 %   +sey*diag(offdiag4,N))/dx^2;

%CC = diag(ones(Nsq,1));


Rhs = zeros(2*Nsq,1);
indx = (Nsq+1)/2; %uses an odd number of points
%indx = 1; %in the corner
%Rhs(indx) = 1;
%Rhs(indx+Nsq) = -1;

Rhs(Nsq+1:2*Nsq) = -1/(Nsq-1); % remove current everywhere else
Rhs(indx+Nsq) = 1; % add current at the center
 
indx1 = Nsq + [1:N];
indx2 = Nsq + (N-1)*N+[1:N];
indx3 = Nsq + [1:N]*N - N + 1;
indx4 = Nsq + [1:N]*N;
cout = 1/(4*N-4);
%Rhs(indx1) = -cout;
%Rhs(indx2) = -cout;
%Rhs(indx3) = -cout;
%Rhs(indx4) = -cout;

% remove current from the boundaries

% change the first row of M
%M(1,:) = [zeros(1,Nsq),ones(1,Nsq)];
%Rhs(1) = 0;

%ans1 = M\Rhs;
%ans1 = gmres(M,Rhs,10,1e-6,[],[],[],ans1);

% use Matlab's conjugate gradient squared method  cgs to solve the equations.
ans1 = cgs(M,Rhs,1e-6,3*N);

phi_i = ans1(1:Nsq);
phi_e = ans1(Nsq+1:2*Nsq);

phi_i = reshape(phi_i,N,N);
phi_e = reshape(phi_e,N,N);
phi = phi_i-phi_e;

figure(1)
    K = 16;
    vmax = max(max(phi));
    vmin = min(min(phi));
    rm = min(abs(vmax),abs(vmin));
    dv = 4*rm/K;
    v = -2*rm +[0:K]*dv;
    contour(phi,v,'linewidth',1.5)
    colorbar
    text(50,52,'+','fontsize',18)
    text(50,63,'-','fontsize',20)
    text(50,39,'-','fontsize',20)
    box off

    % save for external plotting
    contoursave = contour(phi,v);
    writematrix(contoursave','contours.dat')

figure(2)
    contour(phi_i,20,'linewidth',1.5);
    colorbar

figure(3) 
    contour(phi_e,20,'linewidth',1.5)
    colorbar



