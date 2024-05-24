
function exit_time
clear all
global dG tt

Ntot=100; dGlow=0; dGhigh=5;
deltaG = linspace(dGlow,dGhigh,Ntot);
tt=1;

for count=1:Ntot
    dG=deltaG(count);

    tau_off = 0; tau_on=0; Ns=20;
    for i=1:Ns
           s=(i-1)/(Ns-1);
           tau_off = tau_off + (1/tt)*exp(UU(s)/tt)*inner_int(s)*(1/Ns);
    end
    koff(count)=1/tau_off;
    koff_approx(count) = ((-UUpp(0)*UUpp(1))^0.5/(tt*pi))*exp(-dG/tt);

    for i=1:Ns
           s=(i-1)/(Ns-1);
           tau_on = tau_on + (1/tt)*exp(UU_on(s)/tt)*inner_int_on(s)*(1/Ns);
    end
    kon(count)=1/tau_on;
    kon_approx(count) =  ((-UUpp_on(0)*UUpp_on(1))^0.5/(tt*pi))*exp(-2*dG/tt);
end

figure(1)
    subplot(2,2,1); semilogy(deltaG,koff,deltaG,koff_approx,'red')
    xlabel('Delta G')
    ylabel('koff')
    subplot(2,2,2); semilogy(deltaG,kon,deltaG,kon_approx,'red')
    xlabel('Delta G')
    ylabel('kon')
    fiddle=((-UUpp_on(0)*UUpp_on(1))^0.5)/((-UUpp(0)*UUpp(1))^0.5);
    subplot(2,2,3); semilogy(deltaG,kon./koff,deltaG,fiddle*exp(-deltaG),'red')
    xlabel('Delta G')
    ylabel('kon/koff')


dump = [deltaG' koff' koff_approx' kon' kon_approx' kon'./koff' fiddle*exp(-deltaG)'];
save('delG.dat','-ascii','dump')


Ntot=100; templow=0.03; temphigh=10;
dG = 1;
temp=linspace(templow,temphigh,Ntot);

for count=1:Ntot
    tt=temp(count);

    tau_off = 0; tau_on=0; Ns=20;
    for i=1:Ns
           s=(i-1)/(Ns-1);
           tau_off = tau_off + (1/tt)*exp(UU(s)/tt)*inner_int(s)*(1/Ns);
    end
    koff(count)=1/tau_off;
    koff_approx(count) = ((-UUpp(0)*UUpp(1))^0.5/(pi))*exp(-dG/tt);

    for i=1:Ns
           s=(i-1)/(Ns-1);
           tau_on = tau_on + (1/tt)*exp(UU_on(s)/tt)*inner_int_on(s)*(1/Ns);
    end
    kon(count)=1/tau_on;
    kon_approx(count) =  ((-UUpp_on(0)*UUpp_on(1))^0.5/(pi))*exp(-2*dG/tt);
end

figure(2)
subplot(1,2,1); semilogy(1./temp,koff,1./temp,koff_approx,'red')
xlabel('1/T')
ylabel('koff')
subplot(1,2,2); semilogy(1./temp,kon,1./temp,kon_approx,'red')
xlabel('1/T')
ylabel('kon')

%dump = [temp' koff' koff_approx' kon' kon_approx'];
%save('delT.dat','-ascii','dump')

%%%%%%%%%%%%%%%%%%%%%%%
function out=inner_int(s)
global dG tt
out=quad(@(x)exp(-UU(x)/tt),-10,s);
%%%%%%%%%%%%%%%%%%%%%%%
function out=inner_int_on(s)
global dG tt
out=quad(@(x)exp(-UU_on(x)/tt),-10,s);
%%%%%%%%%%%%%%%%%%%%%%
function out=UU(x)
global dG tt
out = dG*( 0.1319*(x.^6) - 0.0417*(x.^5) - 0.5347*(x.^4) - 1.3333*(x.^3) + 2.7778*(x.^2) );
%%%%%%%%%%%%%%%%%%%%%%
function out=UUpp(x)
global dG tt
out = dG*( 6*5*0.1319*(x.^4) - 5*4*0.0417*(x.^3) - 4*3*0.5347*(x.^2) - 3*2*1.3333*(x.^1) + 2*2.7778 );
%%%%%%%%%%%%%%%%%%%%%%
function out=UU_on(x)
global dG tt
out = UU(-x+2)+dG;
%%%%%%%%%%%%%%%%%%%%%%
function out=UUpp_on(x)
global dG tt
out = UUpp(-x+2);
%%%%%%%%%%%%%%%%%%%%%%
