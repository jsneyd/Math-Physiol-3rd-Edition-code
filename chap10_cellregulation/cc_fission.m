clear all

global k1 k2p k2pp k2ppp k3p k3pp k4p k4 k5p k5pp k6 k7 k8
global J3 J4 J5 J7 J8 k9 k10 mu Vawee Viwee Jawee Jiwee 
global Va25 Vi25 Ja25 Ji25 kweep kweepp k25p k25pp
global k11 k12 k13 k14 k15 k16 k12pp Kdiss k16p k16pp
global J9 J10 J11 J12 J13 J14 J15 J16

k1=0.03;k2p=0.03;k2pp=1.0;k2ppp=0.1
k3p=1.0;k3pp=10.0;J3=0.01;k4p=2.0;k4=35.0;J4=0.01
k5p=0.005;k5pp=0.3;k6=0.1;J5=0.3
k7=1.0;k8=0.25;J7=0.001;J8=0.001
k9=0.1;k10=0.04;J9=0.01;J10=0.01
k11=0.1;k12=0.01;k12p=1.0;k12pp=3.0;Kdiss=0.001
k13=0.1;k14=0.1
k15=1.5;k16p=1.0;k16pp=2.0;J15=0.01;J16=0.01
Vawee=0.25;Viwee=1.0;Jawee=0.01;Jiwee=0.01
Va25=1.0;Vi25=0.25;Ja25=0.01;Ji25=0.01
kweep=0.15;kweepp=0.6;k25p=0.05;k25pp=5.0
mu=0.005

npts=50; tend=1;
tfinal=200; test=0;

y0(1)=0.0204; y0(2)=0.9148; y0(3)=0.05;
y0(4)=0.0001; y0(5)=0.062; y0(6)=0.4943;

for loop=1:tfinal
    keep(loop,1)=(loop-1)*tend;     % Keep the times here
    keep(loop,2:7)=y0;              % Keep the solutions here
    tspan=linspace(0,tend,npts);
    [T,Y] = ode45(@cc_fissionodes,tspan,y0);
    if (Y(npts,1)>0.4) test=1; end
    if ((Y(npts,1)<0.05) & (test==1)) 
        Y(npts,6)=Y(npts,6)/2.0;  % Reset the mass for cell division
        test=0;
    end
    y0=Y(npts,1:6);
end



plot(keep(:,1),keep(:,2))
hold on;
plot(keep(:,1),keep(:,7))
hold off;
save gen6.dat keep -ascii

