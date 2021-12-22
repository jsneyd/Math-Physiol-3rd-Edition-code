
clear all

w=linspace(0.6,1.2,500);

mu=0;
w0=1;

a=0.1;
U1=mu^2 + (w0-w).^2;
U2=-mu*mu + 3*(w0-w).^2;
D=27*a*a + 16*(mu^3) - 18*mu*U1;
s=D+(D.*D+4*(U2.^3)).^0.5;
s=s./2;
bigA1= s.^(1/3) + 2*mu - U2./(s.^(1/3));
bigA1 = (bigA1./3).^0.5;

a=0.01;
U1=mu^2 + (w0-w).^2;
U2=-mu*mu + 3*(w0-w).^2;
D=27*a*a + 16*(mu^3) - 18*mu*U1;
s=D+(D.*D+4*(U2.^3)).^0.5;
s=s./2;
bigA2= s.^(1/3) + 2*mu - U2./(s.^(1/3));
bigA2 = (bigA2./3).^0.5;

a=0.0001;
U1=mu^2 + (w0-w).^2;
U2=-mu*mu + 3*(w0-w).^2;
D=27*a*a + 16*(mu^3) - 18*mu*U1;
s=D+(D.*D+4*(U2.^3)).^0.5;
s=s./2;
bigA3= s.^(1/3) + 2*mu - U2./(s.^(1/3));
bigA3 = (bigA3./3).^0.5;


plot(w,bigA1)
hold on
plot(w,bigA2)
plot(w,bigA3)
hold off

dum=[w' bigA1' bigA2' bigA3'];
save('chinchilla.dat','-ascii','dum')