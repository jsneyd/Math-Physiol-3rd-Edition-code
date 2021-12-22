clear
global aA gA aB gB aP gP aM gM aL gL bA bL K K1 KA KL KL1 KLe Le

%specify parameter values
aA = 1.76;
gA=0.52;
aB = 1.66E-2;
gB=2.26e-2;
aP=10;
gP=0.65;
aM =9.97;
gM=0.41;
aL = .288;
gL = 2.26e-2;
bA = 2.15;
bL=2.65e3;
K =6000;
K1=2.52e4;
KA = 1.95;
KL = 9.7e-7;
%KL = 9.7e-4;
KL1 = 1.81;
KLe = 0.26;

% steady state solution
Le = .005; % fig 1
Le = .015; % fig 2  three solutions
%Le = .035; % fig 2
Amin = 0.001;
Amax = .1;
A = [0.0001:.001:300];
M1 = aM*(1+K1*A.^2) ./(( K + K1*A.^2)* gM);
tmp = aL *aP*Le/(( KLe + Le)*gP) - bA *aB *A./(( KA + A)* gB);
M2 = gA *A./tmp;


figure(1)

loglog(A,M1,'b',A,M2,'r','linewidth',2)
axis([1e-4 100 1e-4 200])
xlabel('Allolactose','fontsize',16)
ylabel('mRNA','fontsize',16)
text(2,10,'dM/dt = 0','fontsize',16)
%text(.001,.1,'dA/dt = 0','fontsize',16) %fig 1
text(.15,.2,'dA/dt = 0','fontsize',16) %fig 2
