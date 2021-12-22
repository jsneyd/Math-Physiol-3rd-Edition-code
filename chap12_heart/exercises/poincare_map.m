theta = [0:.01:1];
AA = [.5;.9]
for j=1:2
A = AA(j);
alf = 2*pi;

phi = atan(sin(alf*theta)./(A+cos(alf*theta)))/alf +(A+cos(alf*theta)<0)/2 ...
+(A<1).*((A+cos(alf*theta)>0)&(sin(alf*theta)<0));

figure(1)

plot(theta,phi,'linewidth',2)
hold on
end
hold off
text(.6,.62,'A=0.5','fontsize',16)
text(.6,.85,'A=0.9','fontsize',16)
xlabel('\theta/2\pi','fontsize',16)
ylabel('\phi/2\pi','fontsize',16)

AA = [1.5;1.1];
for j=1:2
A = AA(j);
alf = 2*pi;

phi = atan(sin(alf*theta)./(A+cos(alf*theta)))/alf +(A+cos(alf*theta)<0)/2 ...
+(A<1).*((A+cos(alf*theta)>0)&(sin(alf*theta)<0));

figure(2)

plot(theta,phi,'linewidth',2)
hold on
end
hold off
text(.2,.05,'A=1.5','fontsize',16)
text(.2,.15,'A=1.1','fontsize',16)
xlabel('\theta/2\pi','fontsize',16)
ylabel('\phi/2\pi','fontsize',16)

% plot critical curves
y = [0.25:.01:0.75];

A = -cos(alf*y);
T = y -1/4 + (y>=.5)/2;
figure(3)
plot(T,A,'linewidth',2)
text(.05,.8,'1:1','fontsize',16)
text(.9,.8,'1:1','fontsize',16)
text(.45,.8,'no 1:1','fontsize',16)
xlabel('T/2\pi','fontsize',16)
ylabel('A','fontsize',16)

% check this
y = .863;
A = -cos(alf*y);
T = y+1/4;

A = .9;
T = 0.9;
phi =T+ atan(sin(alf*theta)./(A+cos(alf*theta)))/alf +(A+cos(alf*theta)<0)/2 ...
+(A<1).*((A+cos(alf*theta)>0)&(sin(alf*theta)<0));

figure(4)

plot(theta,phi,theta,theta,'--',theta,theta+1,'--','linewidth',2)
