function GL_glycolysis_sim
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

% this code is to simulate the Goldbeter-Lefever enzyme model

%parameters
global nu eta  
 nu=200;
 eta=120;
 

% set up the integration
tspan = [0:0.0005:1];            % time interval for the solution
IC = [37,1];                        % initial condition
[t,U] = ode15s(@(t,y)rhs(t,y),tspan,IC);

figure(1)
plot(t,U(:,1),t,U(:,2),'--','linewidth',2)
legend('\sigma_1','\sigma_2','fontsize',18,'location','northwest')
xlabel('Time','fontsize',18)
ylabel('Concentration','fontsize',18)

figure(2)
plot(U(:,1),U(:,2),'linewidth',2)
xlabel('\sigma_1','fontsize',20)
ylabel('\sigma_2','fontsize',20)
hold on
% add null clines
w=[0.02:0.1:20];


u1=nu./(1+w).^2;

u2 =eta*w./(1+w).^2;
 
plot(u2,w,'--',u1,w,'.','linewidth',2)
axis([0 40 0 20])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=rhs(t,y)
global nu eta 
u=y(1);
w=y(2);
f =u*(1+w)^2;
Fu = nu-f;
Fw = f -eta*w;
 
out = [Fu Fw]';
