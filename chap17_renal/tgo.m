
function main
clear all
global numx delx delt numtimesteps delay 
global c K1 K2 Cop k

numx=100; delx=1/numx; delt=0.001;
numtimesteps=20000; delay=200;
keepend=zeros(numtimesteps+delay,1);

c=zeros(numx,1);
c(1)=1;

gamma=2;
k=1; Cop=exp(-k); 
K1=1.5; K2=gamma/(K1*k*exp(-k));


for i=1:numtimesteps
    keepend(i+delay)=c(numx);
    c(2:numx) = c(2:numx) - ...
                delt*k*c(2:numx) - phi(keepend(i))*(delt/delx)*(c(2:numx)-c(1:numx-1));
end

figure(1)
plot(keepend)
save('tgoosc.dat','-ascii','keepend')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=phi(x)
global c K1 K2 Cop k

out = 1+K1*tanh(K2*(Cop-x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

