close all; clear all; clc;
global SI r1 EG0 R0 sigma alpha k d0 r2
SI=1;
r1=1;
EG0=1.44;
R0=864;
sigma=43.2;
alpha=20000;
k=432;
d0=0.06;
r2=0.24e-5;

%% fast-slow phase plane
G = linspace(50,1000,1000);
beta = (alpha + G.*G)./(sigma*G.*G) .* k/SI .*(R0./G - EG0);
plot(beta,G)
xlim([0,100])
xlabel('\beta')
ylabel('G')
hold on
dlmwrite('topp.dat',[beta' G']);

%% solve the ode and plot solution
tspan=linspace(0,100,10000);
[t,sol] = ode15s(@(t,y)rhs(t,y),tspan,[100,100,12]);
plot(sol(:,3),sol(:,1),'LineWidth',2)
dlmwrite('topp.dat',[sol(:,3) sol(:,1)],'-append','delimiter',' ','roffset',1);
[t,sol] = ode15s(@(t,y)rhs(t,y),tspan,[100,100,17]);
plot(sol(:,3),sol(:,1),'LineWidth',2)
dlmwrite('topp.dat',[sol(:,3) sol(:,1)],'-append','delimiter',' ','roffset',1);
[t,sol] = ode15s(@(t,y)rhs(t,y),tspan,[400,100,50]);
plot(sol(:,3),sol(:,1),'LineWidth',2)
dlmwrite('topp.dat',[sol(:,3) sol(:,1)],'-append','delimiter',' ','roffset',1);
[t,sol] = ode15s(@(t,y)rhs(t,y),tspan,[500,100,25]);
plot(sol(:,3),sol(:,1),'LineWidth',2)
dlmwrite('topp.dat',[sol(:,3) sol(:,1)],'-append','delimiter',' ','roffset',1);
[t,sol] = ode15s(@(t,y)rhs(t,y),tspan,[500,100,15]);
plot(sol(:,3),sol(:,1),'LineWidth',2)
dlmwrite('topp.dat',[sol(:,3) sol(:,1)],'-append','delimiter',' ','roffset',1);
[t,sol] = ode15s(@(t,y)rhs(t,y),tspan,[600,100,10]);
plot(sol(:,3),sol(:,1),'LineWidth',2)
dlmwrite('topp.dat',[sol(:,3) sol(:,1)],'-append','delimiter',' ','roffset',1);

%% steady states

G1 = (r1/1000 + sqrt((r1/1000)^2 - 4*r2*d0))/(2*r2)
G2 = (r1/1000 - sqrt((r1/1000)^2 - 4*r2*d0))/(2*r2)


%% the ODEs
function out = rhs(t,y)
global SI r1 EG0 R0 sigma alpha k d0 r2
G=y(1);
I=y(2);
beta=y(3);
out = [R0 - (EG0+SI*I)*G
       beta*sigma*G*G/(alpha+G*G) - k*I
       (-d0+r1*G/1000-r2*G*G)*beta];
end

 