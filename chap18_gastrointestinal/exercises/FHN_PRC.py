
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
plt.rcParams['axes.spines.right'] = False  # Me being pedantic and removing the box around the graph
plt.rcParams['axes.spines.top'] = False

def deriv(t, y):
    """ODEs for FHN."""
    V,w = y
    alpha = 0.05
    eps = 0.01
    I = 5
    perturbation = 0
    if (t>tstim) and (t<tstim+0.06):
        perturbation = stim
    Vdot = (1/eps)*( V*(1-V)*(V-alpha) - w) + I + perturbation
    wdot = V - w
    return [Vdot,wdot]


# Output times and range
t0=0
tf=20
t = np.arange(t0,tf,0.001)
# Initial conditions
y0 = [0.1,0]

# No stimulus
stim = 0
tstim = 0 # just a dummy variable, as we need some value to solve the ode
soln = solve_ivp(deriv, (t0, tf), y0, method='Radau',t_eval=t)
V,w = soln.y
peaks,_ = find_peaks(V)
tpeaks = t[peaks]
period = (tpeaks[-1] - tpeaks[-2])/1


# With stimulus
n = 150 #number of points for calculating the phase change
phase_change = np.empty(n)
timestim = np.linspace(tpeaks[7],tpeaks[8],n) # pick the places to do the impulse
stim = 1
for i in range(n):
      tstim = timestim[i]
      soln = solve_ivp(deriv, (t0, tf), y0, method='Radau',t_eval=t)
      peaks2,_ = find_peaks(soln.y[0])
      tpeaks2 = t[peaks2] 
      phase_change[i] = (tpeaks[-1] - tpeaks2[-1])/period  # phase change calculated against last peak. Easier.
    
plt.plot(timestim,phase_change)
np.savetxt('FHN_PRC.dat',np.c_[timestim/period,phase_change])
