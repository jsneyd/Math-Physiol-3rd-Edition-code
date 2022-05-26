import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


def deriv(x,n,v):
    ndot = -(1/v)*(f(x)*(1-n) - g(x)*n)   
    return ndot

def f(x):
    return np.piecewise(x,[x<=0,0<=x<=1,x>1],[0,f1*x,0])
    
def g(x):
    return np.piecewise(x,[x<=0,x>0],[g2,g1*x])

def get_load(v):
    soln = solve_ivp(lambda x,n: deriv(x,n,v),[2, -5],[n0],method='LSODA',max_step=0.005)
    return np.trapz(soln.t,soln.t*soln.y[0])

# ----------------------------------------------------------------

f1=43.3
g1=10 
g2=209
vv = [1,10,50,100,120]

# plot the n for a range of v
# Initial conditions
n0 = 0
plt.figure(1)
for i in range(0,5):
    v = vv[i]
    # You have to prevent the solver taking big steps. It wants to. It's very naughty.
    soln = solve_ivp(lambda x,n: deriv(x,n,v),[2, -3],[n0],method='LSODA',max_step=0.005)
    plt.plot(soln.t,soln.y[0])
plt.xlabel('$x')
plt.ylabel('$n(x)$')

# next, plot the force-velocity curve
p = np.empty(20)
vv = np.empty(20)
for i in range(0,20):
    v = 0.1 + i*120/20
    p[i] = get_load(v)
    vv[i] = v

plt.figure(2)
plt.plot(p,vv)
plt.xlabel('$p$')
plt.ylabel('$v$')

# Plot the maximal velocity as a function of g2
gg2 = np.empty(20)
vmax = np.empty(20)
for i in range(0,20):
    g2 = 150 + i*100/20
    gg2[i] = g2
    vmax[i] = fsolve(get_load,100)
plt.figure(3)
plt.plot(gg2,vmax)
plt.xlabel('$g_2$')
plt.ylabel('$v_{max}$')
    
# Plot the maximal velocity as a function of g1
# first reset the values
f1=43.3
g1=10 
g2=209
gg1 = np.empty(20)
vmax = np.empty(20)
for i in range(0,20):
    g1 = 5 + i*10/20
    gg1[i] = g1
    vmax[i] = fsolve(get_load,107)
plt.figure(4)
plt.plot(gg1,vmax)
plt.xlabel('$g_1$')
plt.ylabel('$v_{max}$')
plt.ylim([105,110])

# Plot the maximal velocity as a function of f1
# first reset the values
f1=43.3
g1=10 
g2=209
ff1 = np.empty(20)
vmax = np.empty(20)
for i in range(0,20):
    f1 = 20 + i*80/20
    ff1[i] = f1
    vmax[i] = fsolve(get_load,107)
plt.figure(5)
plt.plot(ff1,vmax)
plt.xlabel('$f_1$')
plt.ylabel('$v_{max}$')
plt.ylim([105,110])
    
