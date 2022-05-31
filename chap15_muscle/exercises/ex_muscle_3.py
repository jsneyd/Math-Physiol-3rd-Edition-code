
from sympy import *

x= symbols('x')
h,lam,fmax,v,alpha= symbols('h lam fmax v alpha',positive='true')
expr = (exp((x-h)/lam) - alpha)*(exp((x-h)/lam) - exp(fmax*(x-h)/v) )
ans = integrate(expr,(x,-oo,h))
ans1 = simplify( lam/((1-alpha)*(fmax*lam-v))*ans  )

print(ans1)
