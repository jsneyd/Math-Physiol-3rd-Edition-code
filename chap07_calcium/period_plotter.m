% make  diagram from XPP period data

load period_0.dat
A = period_0(2:end,:);
ndxa = find(A(:,4)==2);

load period_0p05.dat
B = period_0p05(2:end,:);
ndxb = find(B(:,4)==2);

load period_0p1.dat
C = period_0p1(2:end,:);
ndxc = find(C(:,4)==3);

load period_0p2.dat
D = period_0p2(2:end,:)
ndxd = find(D(:,4)==2);
figure(1)

plot(A(ndxa,1),A(ndxa,2),'--',B(ndxb,1), B(ndxb,2),C(ndxc,1),C(ndxc,2),D(ndxd,1),D(ndxd,2))

axis([0 2 0 90])

legend('K_{PLC}=0','K_{PLC}=0.05','K_{PLC}=0.1','K_{PLC}=0.2')
