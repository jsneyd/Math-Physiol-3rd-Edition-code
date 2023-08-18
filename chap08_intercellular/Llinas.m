 set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);

 k10=2; k20=1; z1=1; s0=100; PCa=0.00001; F=96490; R=8.315; T=300;
	ci=0; ce=40; n=5;
	Vss=-0.07;
	k1=k10*exp(F*z1*Vss/(R*T));
	ohss=s0*k1/(k1+k20);
	
	x = linspace(0,3,100);

	V=-0.02;
	phi=2*F*V/(R*T);
	k1=k10*exp(F*z1*V/(R*T));
	jj=PCa*2*F*phi*(ci-ce*exp(-phi))/(1-exp(-phi));
	oh_1=k1*(s0/(k1+k20))*(1-exp(-(k1+k20)*x))+ohss*exp(-(k1+k20)*x);
	ICa_1=jj*(s0/n)*((oh_1/s0).^n);

	V=-0.01;
	phi=2*F*V/(R*T);
	k1=k10*exp(F*z1*V/(R*T));
	jj=PCa*2*F*phi*(ci-ce*exp(-phi))/(1-exp(-phi));
	oh_2=k1*(s0/(k1+k20))*(1-exp(-(k1+k20)*x))+ohss*exp(-(k1+k20)*x);
	ICa_2=jj*(s0/n)*((oh_2/s0).^n);

	V=0.02;
	phi=2*F*V/(R*T);
	k1=k10*exp(F*z1*V/(R*T));
	jj=PCa*2*F*phi*(ci-ce*exp(-phi))/(1-exp(-phi));
	oh_4=k1*(s0/(k1+k20))*(1-exp(-(k1+k20)*x))+ohss*exp(-(k1+k20)*x);
	ICa_4=jj*(s0/n)*((oh_4/s0).^n);

	V=0.04;
	phi=2*F*V/(R*T);
	k1=k10*exp(F*z1*V/(R*T));
	jj=PCa*2*F*phi*(ci-ce*exp(-phi))/(1-exp(-phi));
	oh_5=k1*(s0/(k1+k20))*(1-exp(-(k1+k20)*x))+ohss*exp(-(k1+k20)*x);
	ICa_5=jj*(s0/n)*((oh_5/s0).^n);


    figure(1)
    plot(x,ICa_1,x,ICa_2,x,ICa_4,x,ICa_5,'LineWidth',2)
    xlabel('time (ms)')
    ylabel('I_{Ca}(pA/(\mum)^2)')
    legend('V=-20mV','V=-10mV','V=20mV','V=40mV')


    x = linspace(-70,70,100);
 	
 	phi_x=2*F*(x/1000)/(R*T);
 	Jss=PCa*2*F*phi_x.*(ci-ce*exp(-phi_x))./(1-exp(-phi_x));
 	Oss=k10*exp(F*z1*(x/1000)/(R*T)).*(1./(k10*exp(F*z1*(x/1000)/(R*T))+k20));
 	ICass=Jss.*(s0/n).*((Oss).^n);

    figure(2)
    yyaxis left
    plot(x,Jss,x,ICass,'LineWidth',2)
    ylabel('I_{Ca}')
    yyaxis right
    plot(x,Oss,'LineWidth',2)
    xlabel('V (mV)')
   % ylabel('\hat{o}')
    ylim([0,1])
    