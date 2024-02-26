% code to examine phase synchrony

N = 30; % size of grid
Nsq=N^2;
Xcent = 12; 
Ycent = 18; % center of the grid
scale = 50;
amp = 0.3;

[X,Y]= meshgrid(1:N,1:N);
rad= (X-Xcent).^2+(Y-Ycent).^2;
P = rand(N,N);
Pf =  2.*P.*(1-P);

smth=exp(-rad/scale)+Pf;
figure(1)
mxsm = max(max(smth));
smth=smth/mxsm;
pcolor(smth)
 

%create the Nsq x Nsq coupling matrix

offdiag = [ones(N-1,1)]; % off-diagonal 
for j = 1:N-1
offdiag = [offdiag;0;ones(N-1,1)];
end
A=diag(offdiag,1)+diag(offdiag,-1)+diag(ones(N*(N-1),1),N)+diag(ones(N*(N-1),1),-N);
Ad = sum(A);
A=A-diag(Ad);

% create the rhs
rhs=reshape(smth,Nsq,1);

%project out the adjoint nullspace
srhs=sum(rhs);
rhs = rhs-ones(Nsq,1)*srhs/Nsq;

% replace last row of A by ones
A(Nsq,:)=ones(1,Nsq);
rhs(Nsq) = 1;

phase=A\rhs;
nphase = -reshape(phase,N,N);
figure(2)
pcolor(nphase)
 




   

