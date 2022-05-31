#%%
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt

plt.rcParams['axes.spines.right'] = False  # Me being pedantic and removing the box around the graph
plt.rcParams['axes.spines.top'] = False

def f(x):
    f = np.piecewise(x,[x<=0,0<=x<=h,x>h],[0,f1*x/h,0])
    return f

def g(x):
    g = np.piecewise(x,[x<=0,x>0],[g2,g1*x/h])
    return g

def deriv(x,y):
    n0,nm1 = y
    n0dot = -(1/v)*(f(x)*(1 - n0 - nm1) - g(x)*n0)
    nm1dot = -(1/v)*(f(x-delta)*(1 - n0 - nm1) - g(x-delta)*nm1)
    return [n0dot,nm1dot]

k=1.0
rho=1.0
h=1 
delta=0.1
f1=0.25 
g1=0.25 
g2=2

# start just above x=h and integrate backwards in x until x=-2, which seems to be far enough.
x0=h+3*delta
xf=-2
x = np.arange(x0,xf,-0.0001)
# Initial conditions
y0 = [0,0]

# Compute and plot solutions for 4 different values of v
v = 0.01
soln1 = solve_ivp(deriv, (x0, xf), y0, method='Radau',t_eval=x)
v = 0.05
soln2 = solve_ivp(deriv, (x0, xf), y0, method='Radau',t_eval=x)
v = 0.1
soln3 = solve_ivp(deriv, (x0, xf), y0, method='Radau',t_eval=x)
v = 0.5
soln4 = solve_ivp(deriv, (x0, xf), y0, method='Radau',t_eval=x)


# Now just plot out the solutions
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,sharex=True,sharey=True)
ax1.plot(soln1.t,soln1.y[0],soln1.t,soln1.y[1])
ax2.plot(soln2.t,soln2.y[0],soln2.t,soln2.y[1])
ax3.plot(soln3.t,soln3.y[0],soln3.t,soln3.y[1])
ax4.plot(soln4.t,soln4.y[0],soln4.t,soln4.y[1])

plt.legend(['$n_0$','$n_{-1}$'],frameon=False)
ax1.text(-1.5, 0.4, '$v=0.01$')
ax2.text(-1.5, 0.4, '$v=0.05$')
ax3.text(-1.5, 0.4, '$v=0.1$')
ax4.text(-1.5, 0.4, '$v=0.5$')

#plt.savefig('exhill1_2_solutions.png',dpi=300)

#%%
# Now calculate the force-velocity curve. Cycle through the values of v, 
# calculating p by an integral at each step

n = 40
force = np.empty(n)
vel = np.empty(n)
for i in range(n):
    v = 0.01+ 2*i/n
    soln = solve_ivp(deriv, (x0, xf), y0, method='Radau',t_eval=x)
    [n0,nm1] = soln.y
    x = soln.t
    force[i] = -trapezoid(x*n0,x) - trapezoid((x-delta)*nm1,x)   # x is in reverse order so you need to change the sign of the integral
    vel[i] = v

plt.plot(vel,force)
plt.xlabel('v')
plt.ylabel('p')
plt.savefig('exhill1_2_force_velocity.png',dpi=300)
