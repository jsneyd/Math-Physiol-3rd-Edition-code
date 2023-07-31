%apd map figure

clear all
close all
clc

%set parameters
A_max = 616;
A = 750;
mu = 170;
DI_min = 100;
BCL = 680;
lim1 = BCL-DI_min;
lim2 = 2*BCL-DI_min;


x1 = [0:1:lim1];
x2 = [lim1:1:lim2];
y1 = BCL-x1;
y2 = 2*BCL - x2;

g1 = A_max - A*exp(-y1/mu);
g2 = A_max - A*exp(-y2/mu);



x = 400;
nit = 70;
%find the alternans fixed point
for j = 1:nit
    
    y = BCL-x;
    g(2*j-1) = A_max - A*exp(-y/mu);
    g(2*j) = g(2*j-1);
    x = g(2*j-1);
end


plot(x1,g1,x2,g2,[x1,x2],[x1,x2],'--',435,435,'k*',607,607,'k*',g(2*nit-6:2*nit-1),g(2*nit-5:2*nit),'k--','linewidth',2)
axis([0 1200 0 1200])
xlabel('APD_n (ms)','fontsize',16)
ylabel('APD_{n+1} (ms)','fontsize',16)
legend('g_1(APD)','g_2(APD)','one to one line',2)

test=[x1' g1'];
save('test.dat','-ascii','test');
test=[g'];
save('test.dat','-ascii','-append','test');
test=[x2' g2'];
save('test.dat','-ascii','-append','test');

x = 400;
nit = 70;
%find the alternans fixed point
for j = 1:nit
    
    y = BCL-x;
    g(2*j-1) = A_max - A*exp(-y/mu);
    g(2*j) = g(2*j-1);
    x = g(2*j-1);
end


    