import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def deriv(t,y):
    dd = np.zeros(4*num)
    x = y[0:num]
    nm = y[num:2*num]  
    nam = y[2*num:3*num]
    namp = y[3*num:4*num]
    nmp = np.ones(num) - nm - nam - namp
    
    dum = np.trapz(x*(ff(x)*nmp - gg(x)*namp - ggL(x)*nam),x)
    v = dum/(1+np.trapz(namp+nam,x))
    dd[0:num] = -v
    dd[num:2*num] = k2*nmp - k1(t)*nm + ggL(x)*nam
    dd[2*num:3*num] = k6*namp - (k5(t)+ggL(x))*nam
    dd[3*num:4*num] = k5(t)*nam + ff(x)*nmp - (k6+gg(x))*namp
    return dd

def f(x):
    return np.piecewise(x,[x<=0,0<=x<=1,x>1],[0,f1*x,0])
    
def g(x):
    return np.piecewise(x,[x<=0,0<=x],[g2,g1*x])

def gL(x):
    return np.piecewise(x,[x<=0,0<=x],[gL2,gL1*x])

def k1(t):
    return 1 + 0.9*np.sin(0.5*t)

def k5(t):
    return 1 + 0.9*np.sin(0.5*t)

# ----------------------------------------------------------------
ff = np.vectorize(f)
gg = np.vectorize(g)
ggL = np.vectorize(gL)

f1  = 0.88; 
g1  = 0.21;
g2  = 4.4;
gL1 = 0.01;
gL2 = 0.2;
h   = 1;
k2  = 0.1;
k6  = 0.1;

num = 50;       # number of space points
numt = 100;     # number of time outputs
tend = 20;      # final time

# Initial conditions
x0 = np.linspace(-2,2,num);
nm0 = np.ones(num)
nam0 = np.zeros(num)
namp0 = np.zeros(num)
y0 = np.concatenate((x0,nm0,nam0,namp0))  # the final ode is for the spring length
tout = np.linspace(0,tend,numt)

# solve the odes
soln = solve_ivp(deriv,[0, tend],y0,method='Radau',t_eval=tout,rtol=10e-6,atol=10e-9)

# calculate and plot the force
plt.figure(1)
p = np.zeros(numt)
for i in range(0,numt):
    x = soln.y[0:num,i]
    nam = soln.y[2*num:3*num,i]
    namp = soln.y[3*num:4*num,i]
    p[i] = np.trapz( x*(nam+namp),x )
plt.plot(soln.t,p)

np.savetxt('test.dat',np.transpose([soln.t,p]),delimiter=' ')
