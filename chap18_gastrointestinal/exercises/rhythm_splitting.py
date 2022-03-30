import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
plt.rcParams['axes.spines.right'] = False  # Me being pedantic and removing the box around the graph
plt.rcParams['axes.spines.top'] = False

def phase_deriv(t, y):
    phi = y 
    phidot = -2-eps-2*np.sin(phi)
    return phidot

eps = 0.01  

tend = 100
t = np.linspace(0,tend,2000)
soln = solve_ivp(phase_deriv, (0, tend), [0], method='Radau',t_eval=t,rtol=1e-8,atol=1e-8)
plt.plot(soln.t,soln.y[0]/np.pi)
plt.xlabel('t')
plt.ylabel('$\phi/\pi$')

  




#
