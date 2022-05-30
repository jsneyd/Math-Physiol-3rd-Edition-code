# Huxley model of muscle pulling against spring. I've used simpler (and continuous) versions
# of f and g because Python really doesn't like the discontinuous versions. Not sure why, and possibly
# that's fixable. May even be because of a bug, but I suspect that the Python ode routines just
# aren't as good as ode15s.

# We don't actually need to include an ode for the spring extension (L), but we do anyway,
# as this provides a useful check on whether or not your integration is accurate.

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def deriv(t,y):
    dd = np.zeros(2*num+1)
    x = y[0:num]
    n = y[num:2*num]   
    v = np.trapz( x*(ff(x)*(1-n) - gg(x)*n)/(1+np.trapz(n,x)),x)
    dd[0:num] = -v
    dd[num:2*num] = ff(x)*(1-n) - gg(x)*n
    dd[2*num] = v
    return dd

def f(x):
    return np.piecewise(x,[x<=0,0<=x<=1,x>1],[0,f1*x,0])
    
def g(x):
    return np.piecewise(x,[x<=0,0<=x],[g2,g1*x])

# ----------------------------------------------------------------
ff = np.vectorize(f)
gg = np.vectorize(g)

f1=43.3
g1=10 
g2=209
num = 250  # number of space points
numt = 150  # number of time outputs
tend = 0.2

# Initial conditions
x0 = np.linspace(-2,2,num);
n0 = np.zeros(num)
L0 = [0]   # initial extension of spring
y0 = np.concatenate((x0,n0,L0))  # the final ode is for the spring length
tout = np.linspace(0,tend,numt)

soln = solve_ivp(deriv,[0, tend],y0,method='BDF',t_eval=tout,rtol=10e-6,atol=10e-9)
plt.figure(1)
plt.plot(soln.t,soln.y[2*num,:])
plt.xlabel('$t$')
plt.ylabel('$L$')

#recalculate v(t) for each time so you can plot it. Not efficient, just easy
plt.figure(2)
vv = np.zeros(numt)
for i in range(0,numt):
    xx = soln.y[0:num,i]
    nn = soln.y[num:2*num,i]
    vv[i] = np.trapz( xx*(ff(xx)*(1-nn) - gg(xx)*nn)/(1+np.trapz(nn,xx)),xx )
    plt.plot(xx,nn)
    
plt.figure(3)
plt.plot(soln.t,vv)

print('initial velocity = ', np.trapz( x0*(ff(x0)*(1-n0) - gg(x0)*n0)/(1+np.trapz(n0,x0)),x0))
print('maxforce = ',np.trapz(x0*ff(x0)/(ff(x0)+gg(x0)),x0))  # the maximum extension
print('final force is ',soln.y[2*num,-1])








    


    
