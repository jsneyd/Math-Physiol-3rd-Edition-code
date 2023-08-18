% facilitation curve

tpcp=200;
k1=3.75e-3;
k2=2.5e-3;
k3=5e-4;

km1=4e-4;
km2=1e-3;
km3=0.1;
K1=km1/k1;
K2 = km2/k2;
K3 = km3/k3;

T=[1/0.186:0.001:10000];
a1=exp(-km1*(T+tpcp/K1));  
a2=exp(-km2*(T+tpcp/K2));
a3=exp(-km3*(T+tpcp/K3));

Fmax = (1./(1-a1)).*(1./(1-a2)).*(1./(1-a3)); 
    

semilogx(1000./T,Fmax)
xlabel('Stimulus Frequency (Hz)')
ylabel('Facilitation')
axis([0.1 100 0 7])