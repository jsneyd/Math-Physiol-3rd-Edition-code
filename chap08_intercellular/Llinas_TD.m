function Llinas
%a time dependent model of calcium current
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 2.0, ...
   'defaultlinelinewidth', 2.0);
global k10 k2 z1 s0 PCa F ci ce n Vss k1  ohss FbyRT Vr
 k10=2; k2=1; z1=1; s0=100; PCa=0.00001; F=96490; R=8.315; Tp=300;
	ci=0; ce=40; n=5;
	Vss=-70;
	FbyRT = F/(R*Tp*1000);  % units of mV^(-1)
    k1=k10*exp(FbyRT*z1*Vss);
	ohss= k1/(k1+k2);
   Vr=-70;
	 
 
	% there are three equations
 

tt=[0:.01:3];
init=[-50, 0,0.11];
[T,sol] = ode15s(@(t,S)rhs(t,S),tt,init);
 V=sol(:,1);
    oh = sol(:,3);
    currents= fCa(V,oh);
    
   ICa=currents;
    figure(1)
    plot(T,V, 'LineWidth',2)
    xlabel('time (ms)')
     ylabel('V (mV)')
 
    yyaxis right
    plot(T,ICa, 'LineWidth',2)
   ylabel('I_{Ca}(pA/(\mum)^2)')
   
   figure(2)
    plot(T,oh)
xlabel('time (ms)')
ylabel('open probability')
    
end
    
    
    
function currents=fCa(V,oh)
global  s0 PCa F ci ce n FbyRT
  
    phi=2*FbyRT*V;
	jj=PCa*2*F.*phi.*(ci-ce.*exp(-phi))./(1-exp(-phi));
	ICa =(s0/n).*jj.*(oh.^n);
 
     currents = [ICa];
end

 

    function out=rhs(t,s)
    global  Vr k10 k2 z1 FbyRT

    v1=s(1);
    w=s(2);
    op=s(3);

    % FHN odes for the presynaptic voltage
    v1p=100*(0.0001*(v1-Vr)*(70-(v1-Vr))*((v1-Vr)-7)-w);
    wp =0.25*(v1-Vr-5*w);
    %Linas ode for open probability 
    k1=k10*exp(FbyRT*z1*v1);
   
    op = k1.*(1-op)-k2*op;
    out  =[v1p,wp,op]';
    end
     