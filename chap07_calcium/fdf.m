%this file examines the "fire-diffuse-fire" model
clear all;
klist=[0;.0001;.0005;.001;.005;.01;.05;.1;.5;1]%;5;10;50;100];
J = length(klist);
dkappa = 0.001;
kappa = 0;
tau = [.01:.05:50];
for j = 1:J
kappa = klist(j);
for n = 1:1000
g(n,:) = exp(-n./tau/4).*sqrt(1./(4*pi*n*tau)).*exp(-kappa*n*tau);
end

G(j,:) = sum(g);
K(j) = kappa;
h(j,:) = g(1,:);
end
figure(1)
plot(tau,G(1,:),'linewidth',2)
xlabel('\eta','fontsize',16)
ylabel('g(\eta)','fontsize',16)
figure(2)
plot(G(1,:),tau,'linewidth',2)
ylabel('\eta','fontsize',16)
xlabel('c^*L/\sigma','fontsize',16)

figure(3)
plot(tau,G,'linewidth',2)
xlabel('\eta','fontsize',16)
ylabel('g_\beta(\eta)','fontsize',16)
figure(5)
plot(G,tau,'linewidth',2)
ylabel('delay','fontsize',16)
xlabel('threshold','fontsize',16)
axis([0 1 0 5])


figure(4)
semilogy(sqrt(klist),max(G'),sqrt(klist),exp(-sqrt(klist)),'--','linewidth',2)
xlabel('\beta','fontsize',16)
ylabel('g_{max}','fontsize',16)
