
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def deriv(t, y):
    """ODEs for FHN."""
    V = y[0:n]
    w = y[n:2*n]
    alpha = 0.4
    eps = np.linspace(0.02,0.08,num=n)
    Vdot = np.empty(n)
    
    wdot = V - alpha
    Vdot = (1/eps)*( V*(1-V)*(V-alpha) - w)
    for i in range(1,n-2):
        Vdot[i] = Vdot[i] + G*(V[i+1]-V[i]) + G*(V[i-1]-V[i])
    Vdot[0] = Vdot[0]   + G*(V[1]-V[0]) 
    Vdot[n-1] = Vdot[n-1]  + G*(V[n-2]-V[n-1])
    
    return np.concatenate((Vdot,wdot))


n = 200   # number of oscillators in chain

# Output times and range
t0=0
tf=60
t = np.arange(t0,tf,0.05)
# Initial conditions for all oscillators: V = 0.1; w = 0.
y0 = 0.1*np.ones(2*n)
y0[n:2*n] = 0

# First solve without coupling
G=0
soln = solve_ivp(deriv, (t0, tf), y0, method='Radau',t_eval=t)

# Find the peaks of each oscillator, and thus calculate the frequency
freq = np.empty(n)
for i in range(0, n):
    peaks,_ = find_peaks(soln.y[i])
    peak_times = t[peaks]
    freq[i] = 1/((peak_times[-1] - peak_times[-10])/10)

# Plot solution
plt.plot(freq)

# Now solve with coupling
G=50
soln = solve_ivp(deriv, (t0, tf), y0, method='Radau',t_eval=t)

# Find the peaks of each oscillator, and thus calculate the frequency
freq = np.empty(n)
for i in range(0, n):
    peaks,_ = find_peaks(soln.y[i])
    peak_times = t[peaks]
    freq[i] = 1/((peak_times[-1] - peak_times[-10])/10)

# Plot solution
plt.plot(freq)
