import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def deriv(t,y):
     m0,m1,m2 = y
     a = m0**2/(2*np.sqrt(3)*np.sqrt(m0*m2-m1*m1)) # value of n where non-zero
     dum0 = m1/m0                                  
     dum1 = (np.sqrt(3)/m0)*(np.sqrt(m0*m2-m1*m1))
     a0 = dum0 - dum1                      #  n is non-zero between b0 and b1
     a1 = dum0 + dum1  
     # Now calculate all the limits for the integrals
     (i1c1,i1c2) = get_limits(0,1,a0,a1)
     (i2c1,i2c2) = get_limits(-np.inf,0,a0,a1)
     (i3c1,i3c2) = get_limits(0,np.inf,a0,a1)
     
     (int1,int2,int3) = get_integrals(0,a,i1c1,i1c2,i2c1,i2c2,i3c1,i3c2)
     m0dot = b0 - int1 - int2 -int3
     int1,int2,int3 = get_integrals(1,a,i1c1,i1c2,i2c1,i2c2,i3c1,i3c2)
     m1dot = b1 - v(t)*m0 - int1 - int2 -int3
     int1,int2,int3 = get_integrals(2,a,i1c1,i1c2,i2c1,i2c2,i3c1,i3c2)
     m2dot = b2 -2*v(t)*m1 - int1 - int2 -int3
     return [m0dot,m1dot,m2dot]

def get_limits(a1,b1,a2,b2):  # calculates the overlap of two regions
    c1 = max(a1,a2)
    c2 = min(b1,b2)
    if (c2<c1):
        c1=0
        c2=0
    return [c1,c2]

def get_integrals(lam,a,i1c1,i1c2,i2c1,i2c2,i3c1,i3c2):
     # we have to calculate 3 integrals for each ode. The integral is calculated exactly.    
     i1 = f1*a/(lam+2)*(i1c2**(lam+2) - i1c1**(lam+2))
     i2 = g2*a/(lam+1)*(i2c2**(lam+1) - i2c1**(lam+1))
     i3 = g1*a/(lam+2)*(i3c2**(lam+2) - i3c1**(lam+2))
     return [i1,i2,i3]

def f(x):
    return np.piecewise(x,[x<=0,0<=x<=1,x>1],[0,f1*x,0])
    
def g(x):
    return np.piecewise(x,[x<=0,x>0],[g2,g1*x])

def v(t):
    return v0 + v1*np.sin(50*t)

# ----------------------------------------------------------------

f1=43.3
g1=10 
g2=209
v0 = 0 # +/- 10 for constant velocity solutions
v1 = 25  # 25 for oscillatory solutions

ff = np.vectorize(f)
gg = np.vectorize(g)

# This is the range of x for which we do all the integrations
xp = np.linspace(-1,3,num = 1000)
b0 = np.trapz(ff(xp),x=xp)
b1 = np.trapz(xp*ff(xp),x=xp)
b2 = np.trapz(xp*xp*ff(xp),x=xp)


# Initial conditions
n0 = ff(xp)/(ff(xp)+gg(xp))
m0_0 = np.trapz(n0,x=xp)
m1_0 = np.trapz(xp*n0,x=xp)
m2_0 = np.trapz(xp*xp*n0,x=xp)

tp = np.linspace(0,0.5,100)
soln1 = solve_ivp(deriv,[0, 0.5],[m0_0,m1_0,m2_0],method='Radau',t_eval=tp)
plt.plot(soln1.t,soln1.y[0],soln1.t,soln1.y[1],soln1.t,soln1.y[2])

#np.savetxt('moments',np.transpose(soln1.y))
#np.savetxt('moments_time',np.transpose(soln1.t))
