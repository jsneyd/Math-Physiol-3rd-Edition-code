clear all

global k1 k2p k2pp k2ppp k3p k3pp k4p k4 k5p k5pp k6 k7 k8
global J3 J4 J5 J7 J8 k9 k10 mu

k1=0.04; k2p=0.04; k2pp=1.0; k2ppp=1.0;
k3p=1.0; k3pp=10.0; k4p=2.0; k4=35.0;
k5p=0.005; k5pp=0.2; k6=0.1; k7=1.0; k8=0.5;
J3=0.04; J4=0.04; J5=0.3;
J7=1.0e-3; J8=1.0e-3;
k9=0.1; k10=0.02;
mu=0.005;

npts=50; tend=1;
tfinal=800; test=0;

y0(1)=0.0204; y0(2)=0.9148; y0(3)=0.05;
y0(4)=0.0001; y0(5)=0.062; y0(6)=0.4943;

for loop=1:tfinal
    keep(loop,1)=(loop-1)*tend;     % Keep the times here
    keep(loop,2:7)=y0;              % Keep the solutions here
    tspan=linspace(0,tend,npts);
    [T,Y] = ode45(@ccgen5odes,tspan,y0);
    if (Y(npts,1)>0.4) test=1; end
    if ((Y(npts,1)<0.05) & (test==1)) 
        Y(npts,6)=Y(npts,6)/2.0;  % Reset the mass for cell division
        test=0;
    end
    y0=Y(npts,1:6);
end



plot(keep(:,1),keep(:,2),keep(:,1),keep(:,7),'linewidth',2)
legend('[Cdc13_T]','m','fontsize',16)
 
save gen6.dat keep -ascii
