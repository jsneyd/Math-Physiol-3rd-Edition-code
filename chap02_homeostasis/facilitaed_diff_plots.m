% facilitated diffusion plots

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

 

%parameters
rho = 10;
s0 = 2;
sL = 0.1;
mu = 1/((1+s0)*(1+sL));

sig = [0:0.1:20];

u =  sig./(1+sig);

sig0 = s0+rho*s0/(1+s0);
fac = (1+mu*rho)*(s0-sL);
xi = -(sig+rho*u-sig0)/fac;
stot = sig+u;
sigp = fac*(1+sig).^2./(sig.^2+rho+2*sig+1);
up = sigp./(1+sig).^2;


figure(1)
plot(xi,sig,'r',xi,u,'b',xi,stot,'g','linewidth',2)
axis([0 1 0 3])
xlabel('y','fontsize',20)
ylabel('Oxygen Concentration','fontsize',20)
text(0.1,2.4,'Total','fontsize',20)
text(0.1,0.45,'Bound','fontsize',20)
text(0.1,1.3,'Free','fontsize',20)
if(cmykflg==1)
    print('../../figs_c/chapt_6/facil_conc_fig','-deps','-cmyk')
end
figure(2)
plot(xi,sigp,xi,rho*up,'--','linewidth',2)
text(0.3,5,'Bound Oxygen Flux','fontsize',20)
text(0.3,2.6,'Free Oxygen Flux','fontsize',20)
axis([0 1 0 7])
xlabel('y','fontsize',20)
ylabel('Oxygen Flux','fontsize',20)
if(cmykflg==1)
    print('../../figs_c/chapt_6/facil_flux_fig','-deps','-cmyk')
end
rho = 5;
gam = sig+rho*u;

figure(3)
plot(gam,sig,sig,sig,'--','linewidth',2)
axis([0 10 0 10])
xlabel('Oxygen Consumption','fontsize',20)
ylabel('Critical external oxygen concentration','fontsize',20)
text(6,2.1,'\rho = 5','fontsize',20)
text(6,5.8,'\rho = 0','fontsize',20)
if(cmykflg==1)
    print('../../figs_c/chapt_6/critcl_oxy_mscle','-deps','-cmyk')
end
r = rho;
g = 14;
xi = sqrt(gam/g);
xi0 = sqrt(sig/g);
s1 = g;
s0 = (-r-s1-1+sqrt(4*r*s1^2+r^2+6*r*s1+s1^2+2*r+2*s1+1))/(1+s1)/2;
sig0=[s0:0.1:20];
gam0 = sig0+rho*sig0./(1+sig0);

xi1 = sqrt(abs(gam0-s0-r*s0/(1+s0))/g);
figure(4)

plot(xi1,sig0,xi,sig,'g',xi0,sig,'r--','linewidth',2)
xlabel('Radius','fontsize',20)
ylabel('Free Oxygen','fontsize',20)
legend('\rho=5','\rho=5, critical external concentration','\rho=0, critical external concentration', ...
    'Location','northwest')
,
axis([0 1 0 14])

if(cmykflg==1)
    print('../../figs_c/chapt_6/free_oxy_mscle','-deps','-cmyk')
end
