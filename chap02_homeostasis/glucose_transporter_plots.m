%glucose_transporter_plots

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
 
k = 0.5;
 se=[0:.1:10];
si=0;
J0 = flux(se,si,k)
si=1;
J1 = flux(se,si,k)
si=2;
J2 = flux(se,si,k)
figure(1)
plot(se,J0,'g',se,J1,'r',se,J2,'b',se,zeros(length(se),1),'--','linewidth',2)
legend('\sigma_i=0','\sigma_i=1','\sigma_i=2','fontsize',18)
xlabel('External glucose','fontsize',20)
ylabel('Glucose flux','fontsize',20)


function J=flux(se,si,k)
 J=(se-si)./((si+1+k).*(se+1+k)-k^2);

end