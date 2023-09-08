function Selkov_glycolysis_sim
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
 
% this code is to simulate the Selkov enzyme model

%parameters
 
par.nu=0.0285;
 

par.eta=0.08;
 

par.alph=1.0;
par.g=2;

% set up the integration
tspan = [0:0.1:800];            % time interval for the solution
IC = [0.4,0.4];                        % initial condition
[t,U] = ode15s(@(t,y)rhs(t,y,par),tspan,IC);

% plot the solution as a function of time:
figure(1)
plot(t,U(:,1),'r',t,U(:,2),'b--','linewidth',2)
legend('boxoff')
legend('\sigma_1','\sigma_2','fontsize',18,'location','northwest')
xlabel('time','fontsize',18)
ylabel('concentration','fontsize',18)

% plot the phase portrait
figure(2)

% add null clines
s2=[0.02:0.01:1.4];
s1a=par.nu/(1-par.nu)*(1+s2.^par.g)./s2.^par.g;
s1b=(1+s2.^par.g)./(s2.^(par.g-1).*(par.alph/par.eta-s2));
plot(U(:,1),U(:,2),'r',s1a,s2,'b--',s1b,s2,'g--','linewidth',2)
legend('boxoff')
legend('solution','\sigma_1 nullcline','\sigma_2 nullcline','fontsize',20)
xlabel('\sigma_1','fontsize',20)
ylabel('\sigma_2','fontsize',20)
axis([0 1.4 0 1])
 eta = ones(51,1)*[0:.01:1];
nu = [0:.01:0.5]'*ones(1,101);

F=2*eta.^3.*nu + nu.^4 - eta.^3 + eta.*nu.^2 - 2*nu.^3 + nu.^2;
figure(3)
contour(nu,eta,F,[0 0],'linewidth',2)
hold on
plot(par.nu,par.eta,'*')
xlabel('\nu')
ylabel('\eta')
text(0.1,0.6,'unstable','fontsize',20)
text(0.3,0.2,'stable','fontsize',20)
hold off

%out = [t U];
%save('test1.dat','out')
%out = [s2' s1a' s1b'];
%save('test2.dat','out')

%>>>>>>> 53bbebdd77194aebd9a7a48e35fd001590167bd9

function out=rhs(t,y,par)
 
u=y(1);
w=y(2);
g=par.g;
f =u*w^g/(w^g*u+w^g+1);
Fu = par.nu-f;
Fw = par.alph*f -par.eta*w;
 
out = [Fu Fw]';
