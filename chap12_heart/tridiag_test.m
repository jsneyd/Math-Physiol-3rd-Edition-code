% this is  to build and test a tridiagonal solver

% The matrix
N = 10;
D0 = -ones(N,1);
Dm1 = .4*ones(N-1,1);
D1 = .3*ones(N-1,1);

A= diag(D0) + diag(Dm1,-1) + diag(D1,1);

Rhs = ones(N,1);

x = A\Rhs;
plot(x)

% now the tridiagonal bit
%forward sweep
m(1) = 1;
entry(1) = D0(1);
for j = 2:N
    m(j) = Dm1(j-1)/entry(j-1);
    entry(j) = D0(j) - m(j)*D1(j-1);
    
end
u(1) = Rhs(1);
% now the forward sweep
for j = 2:N
    u(j) = Rhs(j)-m(j)*u(j-1);
end

% now the backward sweep
u(N) = u(N)/entry(N);
for j = N-1:-1:1
    
    u(j) = (u(j)-D1(j)*u(j+1))/entry(j);
end

plot([1:N],x,'r',[1:N],u,'b')
figure(1)

