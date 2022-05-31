# solve the Huxley model by the method of characteristics for an oscillatory velocity
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def deriv(t,y):
    dd = np.zeros(2*num)
    x = y[0:num]
    n = y[num:2*num]
    dd[0:num] = -v(t)
    dd[num:2*num] = ff(x)*(1-n) - gg(x)*n
    return dd

def f(x):
    return np.piecewise(x,[x<=0,0<=x<=1,x>1],[0,f1*x,0])
    
def g(x):
    return np.piecewise(x,[x<=0,x>0],[g2,g1*x])

def v(t):
    return 50*np.sin(50*t)

# ----------------------------------------------------------------
ff = np.vectorize(f)
gg = np.vectorize(g)

f1=43.3
g1=10 
g2=209
num = 200  # number of space points

# Initial conditions
x0 = np.linspace(-3,4,num);  # Here's the spatial discretisation
n0 = ff(x0)/(ff(x0) + gg(x0))
y0 = np.concatenate((x0,n0))
tend = 0.4
tout = np.linspace(0,tend,11)   # output at these times, for plotting purposes

soln = solve_ivp(deriv,[0, tend],y0,method='LSODA',max_step=0.005,t_eval=tout)
plt.plot(soln.y[0:num,10],soln.y[num:2*num,10])
plt.xlabel('$x$')
plt.ylabel('$n(x)$')

# I've only bothered to plot the output at the final time, but it's easy to plot 
# the outputs at intermediate times. I'll leave that to you.






    


    
