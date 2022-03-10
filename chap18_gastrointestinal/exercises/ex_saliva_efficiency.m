
close all
clear all
clc

a1 = 0.01;
a2 =0.05;

j1=linspace(0,2,100);
j2 = j1;

c1 = a1;
c2 = 1;
c3 = -j1;
q1 = (-c2 + (c2^2 - 4*c1*c3).^0.5)/(2*c1);

c1 = a2;
c2 = 1;
c3 = -j2;
q2 = (-c2 + (c2^2 - 4*c1*c3).^0.5)/(2*c1);

c1 = 1/(1/a1 + 1/a2);
c2 = 1;
c3 = -j1-j2;
Q = (-c2 + (c2^2 - 4*c1*c3).^0.5)/(2*c1);

plot(j1+j2,q1+q2,'r',j1+j2,Q,'b')

%% 

close all
clear all
clc

syms a1 a2 aT j1 j2 jT c1 c2 c3 qq1 qq2 QQ
n=2;

c1 = a1;
c2 = 1;
c3 = -j1;
q1 = (-c2 + (c2^2 - 4*c1*c3).^0.5)/(2*c1);
qq1 = taylor(q1,a1,'Order',n);

c1 = a2;
c2 = 1;
c3 = -j2;
q2 = (-c2 + (c2^2 - 4*c1*c3).^0.5)/(2*c1);
qq2 = taylor(q2,a2,'Order',n);

c1 = aT;
c2 = 1;
c3 = -j1-j2;
Q = (-c2 + (c2^2 - 4*c1*c3).^0.5)/(2*c1);
QQ = taylor(Q,aT,'Order',n);
QQ = subs(QQ,aT,1/(1/a1+1/a2));

pretty(simplify(QQ-qq1-qq2))
%collect(QQ-qq1-qq2,[a1,a2])