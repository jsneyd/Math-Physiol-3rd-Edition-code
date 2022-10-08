function Selkov_glycolysis_sim
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

% this code is to simulate the Selkov enzyme model

%parameters
global nu eta alph g
nu=0.0285
eta=0.1
alph=1.0
g=2

% set up the integration
tspan = [0:0.1:800];            % time interval for the solution
IC = [0.4,0.4];                        % initial condition
[t,U] = ode15s(@(t,y)rhs(t,y),tspan,IC);

% plot the solution as a function of time:
figure(1)
plot(t,U(:,1),t,U(:,2),'--','linewidth',2)
legend('\sigma_1','\sigma_2','fontsize',18,'location','northwest')
xlabel('Time','fontsize',18)
ylabel('Concentration','fontsize',18)

% plot the phase portrait
figure(2)
plot(U(:,1),U(:,2),'linewidth',2)
xlabel('\sigma_1','fontsize',20)
ylabel('\sigma_2','fontsize',20)
hold on
% add null clines
s2=[0.02:0.01:1.4];
s1a=nu/(1-nu)*(1+s2.^g)./s2.^g;
s1b=(1+s2.^g)./(s2.^(g-1).*(alph/eta-s2));
plot(s1a,s2,'--',s1b,s2,'.','linewidth',2)
axis([0 1.4 0 1])
hold off

%out = [t U];
%save('test1.dat','out')
%out = [s2' s1a' s1b'];
%save('test2.dat','out')


function out=rhs(t,y)
global nu eta alph g
u=y(1);
w=y(2);
f =u*w^g/(w^g*u+w^g+1);
Fu = nu-f;
Fw = alph*f -eta*w;
 
out = [Fu Fw]';
