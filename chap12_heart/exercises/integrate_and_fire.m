%This is the integrate and fire model

s0 = 2;
sm = 1;
g = .1;
vt = 10;


a = 2*pi;
t = [0:.01:5];
f = s0/g + sm*(g*sin(a*t) - cos(a*t))/(g^2+1);

g1 = exp(a*g*t).*f;
g2 = exp(a*g*t).*(f-vt);

plot(t,g1,t,g2)
