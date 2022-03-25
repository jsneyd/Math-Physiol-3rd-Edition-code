import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
plt.rcParams['axes.spines.right'] = False  # Me being pedantic and removing the box around the graph
plt.rcParams['axes.spines.top'] = False


# We solve the adjoint equation at the same time as we solve the FHN equations. This way,
# the solutions of the adjoint and of the limit cycle are always lined up in time.
# However, this means you have to be careful to reverse the sign of the derivs
# in the adjoint equations, as they have to be solved backwards in time, for
# stability

# THIS APPROACH DOESN'T WORK

def FHN(t, y):
    v,w = y
    vdot = (1/eps)*( v*(1-v)*(v-alpha) - w) + I
    wdot = v - w
    return [vdot,wdot]

def FHN_adj(t, y):
    Va,Wa = y
    Vadot = -(1/eps)*(-3*v*v + 2*v*(1+alpha) - alpha)*Va - Wa
    Wadot = (1/eps)*Va + Wa
    return [Vadot,Wadot]

alpha = 0.05
eps = 0.01  
I = 5
t = np.linspace(0,20,20000)

# make sure the initial condition is normalised. I don't really know how
# important this is, but it does make the output more consistent, depending
# less on the exact initial condition
v0 = 0.1
w0 = 0
Wa0 = 0
Va0 = (1 - Wa0*(v0-w0))/((1/eps)*(v0*(1-v0)*(v0-alpha) - w0) + I)  
y0 = [v0,w0,Va0,Wa0]

soln = solve_ivp(FHN_adj, (0, 20), y0, method='Radau',t_eval=t)
peaks,_ = find_peaks(soln.y[0])                    
T = t[peaks[-1]] - t[peaks[-2]]                     # calculate the period
# extract two periods of the solution, for plotting and integration
tt = t[peaks[-3]:peaks[-1]]                  
v = soln.y[0,peaks[-3]:peaks[-1]]
w = soln.y[1,peaks[-3]:peaks[-1]]
Va = soln.y[2,peaks[-3]:peaks[-1]]
Wa = soln.y[3,peaks[-3]:peaks[-1]]
Va = np.flip(Va)                        # to get time going in the correct direction
Wa = np.flip(Wa)

plt.figure(1)
plt.plot(tt/T,Va)
#plt.xlim([39,40])
plt.show(1)
np.savetxt('FHN_adjoint.dat',np.c_[tt/T,Va])

# Finally, we have to normalise the solution to the adjoint ODE
# This isn't working yet, as the theory says that the integrand should
# be constant. It clearly isn't.

integrand = Va*((v*(1-v)*(v-alpha) - w)*(1/eps) + I) + Wa*(v-w)
norm = np.trapz(integrand,tt)/(2*T)
#plt.figure(2)
plt.plot(tt/T,integrand)
#plt.xlim([39,40])
plt.show(2)
print(norm)
