
from sympy import *

x, t= symbols('x t',positive = 'true')
rho = 1
k = 1
fplusg = 0.5
g1 = 3/32
g2 = 2
delx = 0.3
h = 1
alpha = (fplusg - g1)/fplusg

dum = x*exp(-fplusg*x*t/h)
ans = integrate(dum,(x,h-delx,h))
print(ans)   # Just for funsies, check the integral

plot((rho*k*alpha*(h**2)/2,(t,-2,0)),
        (rho*k*alpha*(h**2)/2 - alpha*exp(-g2*t)*(delx**2) - alpha*ans,(t,0,5)),
        ylim=(0,0.43),
        axis_center=(0,0),
        ylabel='force',
        xlabel='time')
