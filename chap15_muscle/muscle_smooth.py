
#   -------------------------------------------------------------------
# 
#  Hai-Murphy-Huxley model of smooth muscle. Use the method of
#  characteristics to compute the distributions as functions of space and
#  time. Use two different versions of the model.
# 
#    For Chapter 15, Section 15.8 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

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
    
    dd[0:num] = -v(t)
    dd[num:2*num] = k2*nmp - k1*nm + ggL(x)*nam
    dd[2*num:3*num] = k6*namp - (k5+ggL(x))*nam
    dd[3*num:4*num] = k5*nam + ff(x)*nmp - (k6+gg(x))*namp
    return dd

def deriv_complex(t,y):
    dd = np.zeros(4*num)
    x = y[0:num]
    nm = y[num:2*num]  
    nam = y[2*num:3*num]
    namp = y[3*num:4*num]
    
    # Now do all the shifts of Nam. By hand is the easiest to follow, although it's wordy
    s=xshift # to save typing
    nams1 = np.zeros(num);  nams1[s:-1] = namp[0:-1-s]
    nams2 = np.zeros(num);  nams2[2*s:-1] = namp[0:-1-2*s];
    nams3 = np.zeros(num);  nams3[3*s:-1] = namp[0:-1-3*s];
    namsm1 = np.zeros(num);  namsm1[0:-1-s] = namp[s:-1];  
    namsm2 = np.zeros(num);  namsm2[0:-1-2*s] = namp[2*s:-1];
    namsm3 = np.zeros(num);  namsm3[0:-1-3*s] = namp[3*s:-1];
    nmp = np.ones(num) - nm - namp - nam - nams1 - nams2 - nams3 - namsm1 - namsm2 - namsm3;
    
    dd[0:num] = -v(t)
    dd[num:2*num] = k2*nmp - k1*nm + ggL(x)*nam + (
            ggL(x+delx)*nams1 + ggL(x+2*delx)*nams2 + ggL(x+3*delx)*nams3 +
            ggL(x-delx)*namsm1 + ggL(x-2*delx)*namsm2 + ggL(x-3*delx)*namsm3)
    dd[2*num:3*num] = k6*namp - (k5+ggL(x))*nam
    dd[3*num:4*num] = k5*nam + ff(x)*nmp - (k6+gg(x))*namp
    return dd

def f(x):
    return np.piecewise(x,[x<=0,0<=x<=1,x>1],[0,f1*x,0])
    
def g(x):
    return np.piecewise(x,[x<=0,0<=x],[g2,g1*x])

def gL(x):
    return np.piecewise(x,[x<=0,0<=x],[gL2,gL1*x])

def v(t):
    return np.sin(2*t)

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
delx = 1.2;
k1  = 0.35;
k2  = 0.1;
k5  = 0.35;
k6  = 0.1;

num = 50;      # number of space points
numt = 100;     # number of time outputs
tend = 5;      # final time

# calculate the shifts needed in the complicated model version
xlow = -6; xhigh = 6;
point_spacing = (xhigh-xlow)/num;
xshift = round(delx/point_spacing);        # calculate how many points fit in length delx

# Initial conditions
x0 = np.linspace(xlow,xhigh,num);
nm0 = np.zeros(num)
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

# solve the more complex odes
soln = solve_ivp(deriv_complex,[0, tend],y0,method='Radau',t_eval=tout,rtol=10e-6,atol=10e-9)

# calculate and plot the force
plt.figure(1)
p = np.zeros(numt)
for i in range(0,numt):
    x = soln.y[0:num,i]
    nam = soln.y[2*num:3*num,i]
    namp = soln.y[3*num:4*num,i]
    p[i] = np.trapz( x*(nam+namp),x )
plt.plot(soln.t,p)
