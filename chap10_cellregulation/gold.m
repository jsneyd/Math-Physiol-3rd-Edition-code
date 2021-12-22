function out=gold(a,b,c,d)

out=2.0*a*d/(b-a + b*c + a*d + sqrt((b-a+b*c+a*d)^2 - 4.0*a*d*(b-a)));
