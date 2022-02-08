
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

D = 1000
r = 0.05
c0 = 0.3
P = 0.2
L = 100
a = 0.1
N0 = 0.3

def fun(x, y, p):
    N = N0*(1-np.heaviside(x-a*L,1))
    return np.vstack( (y[1], (1/D)*(2*P*y[0]*(y[0]-c0)/r + y[1]*y[2] - 2*N/r),2*P*(y[0]-c0)/r) )

def bc(ya, yb, p):
    cc = p[0]
    return np.array([ya[1], ya[2], yb[0]-c0,ya[0]-cc])

x = np.linspace(0, 100, num=1000)
y = np.zeros((3, x.size))
y[0] = 1
y[1] = 0
y[2] = 10

sol = solve_bvp(fun, bc, x, y,p=[1])

# All the rest is just plotting out the solution and making it look nicer. 
# The hard work is already done.

x_plot = np.linspace(0, 100, num=1000)

fig, ax1 = plt.subplots(dpi=600)
c_plot = sol.sol(x_plot)[0]
v_plot = sol.sol(x_plot)[2]
ax1.plot(x_plot, c_plot,color='red')
ax1.set_xlabel('x ($\mu$m)')
ax1.set_ylabel('c ($\mu$M)', color = 'red')
ax1.tick_params(axis ='y', labelcolor = 'red')

ax2 = ax1.twinx()
ax2.plot(x_plot, v_plot,color='blue')
ax2.set_xlabel('x ($\mu$m)')
ax2.set_ylabel('v ($\mu$m/s)', color = 'blue')
ax2.tick_params(axis ='y', labelcolor = 'blue')

plt.show()
