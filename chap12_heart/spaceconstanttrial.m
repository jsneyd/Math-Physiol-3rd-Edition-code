%trial space constant code

N=100;
dx=1/N;
A=-diag(2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
A(1,2) = 2;
A(N,N-1) = 2;

R = zeros(N,1);
R(1) = 1;
R(N) = -1;
A=A*N-diag(2*ones(N,1));
vs=A\R;
figure(1)
plot(vs)

