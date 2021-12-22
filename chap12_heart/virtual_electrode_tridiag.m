% this is an attempt to find a virtual electrode in a bidomain calculation
% using an iterative method for the matrix solver
clear

N = 5; % number of grid points (must be odd)
Nsq = N^2;
dx = 5/N;

%equal anisotropy
six = .1;
sex = .2;
siy = 1;
sey = 2;

%conductances 
six = 1;
siy = 0.09;
sex = 0.35;
sey = 0.11;

% nominal  values of Roth
six = 2;
siy - .2;
sex = 8;
sey = 2;

% inverse anisotropy
six = 2;
siy = .2;
sex = .2;
sey = 2;

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

offdiag1 = [ones(N-2,1);2]; % off-diagonal -1
offdiag2 = [2;ones(N-2,1)]; % off-diagonal +1
for j = 1:N-1
    offdiag1 = [offdiag1;0;ones(N-2,1);2];
    offdiag2 = [offdiag2;0;2;ones(N-2,1)];
end

offdiag3 = [ones(N*(N-2),1);2*ones(N,1);zeros(N,1)]; % -N (extended length)
offdiag4 = [2*ones(N,1);ones(N*(N-2),1)]; % +N

%the I matrix
D0i = -2*(six+siy)*ones(Nsq,1);
Dm1i = six*offdiag1;
D1i = six*offdiag2;

Dy1i = siy*offdiag4;
Dym1i = siy*offdiag3;

%the e matrix
D0e = -2*(sex+sey)*ones(Nsq,1);
Dm1e = sex*offdiag1;
D1e = sex*offdiag2;

Dy1e = sey*offdiag4;
Dym1e = offdiag3;

%preprocess the triangular matrices
mi(1) = 1;
me(1) = 1;
entryi(1) = D0i(1);
entrye(1) = D0e(1);

for j = 2:N
    me(j) = Dm1e(j-1)/entrye(j-1);
    entrye(j) = D0e(j) - me(j)*D1e(j-1);
    mi(j) = Dm1i(j-1)/entryi(j-1);
    entryi(j) = D0i(j) - mi(j)*D1i(j-1);
end

%specify the righthand sides

Rhsi = zeros(Nsq,1);
Rhse = zeros(Nsq,1);
indx = (Nsq+1)/2; %uses an odd number of points
%indx = 1; %in the corner
%Rhs(indx) = 1;
Rhse(indx) = -1;
indx1 =  [1:N];
indx2 =  (N-1)*N+[1:N];
indx3 =  [1:N]*N - N + 1;
indx4 =  [1:N]*N;
cout = 1/(4*N-4);
Rhse(indx1) = -cout;
Rhse(indx2) = -cout;
Rhse(indx3) = -cout;
Rhse(indx4) = -cout;

% now set up the iteration
%initial solution
phi_i = zeros(Nsq,1);
phi_e = zeros(Nsq,1);

%the temporary righthand side
for k = 1:10
re = Rhse - phi_i-[Dy1e.*phi_e(N+1:Nsq);zeros(N,1)] - Dym1e.*phi_e;
ri = Rhsi - phi_e - [Dy1i.*phi_i(N+1:Nsq);zeros(N,1)] - Dym1e.*phi_i;

%iterate

% now the forward sweep
for j = 2:N
    re(j) = re(j)-me(j)*re(j-1);
    ri(j) = ri(j)-mi(j)*ri(j-1);
end

% now the backward sweep
re(N) = re(N)/entrye(N);
ri(N) = ri(N)/entryi(N);
for j = N-1:-1:1
    
    re(j) = (re(j)-D1e(j)*re(j+1))/entrye(j);
    
    ri(j) = (re(j)-D1e(j)*re(j+1))/entrye(j);
end
erre = phi_e-re;
erri = phi_i-ri;
max(erre)
max(erri)

phi_i = ri;
phi_e = re;
end

phi_i = reshape(phi_i,N,N);
phi_e = reshape(phi_e,N,N);
phi = phi_i-phi_e;

figure(1)
contour(phi,40)

figure(2)
contour(phi_i,20);

figure(3) 
contour(phi_e,20)




