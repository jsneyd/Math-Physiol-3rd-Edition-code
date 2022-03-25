import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
plt.rcParams['axes.spines.right'] = False  # Me being pedantic and removing the box around the graph
plt.rcParams['axes.spines.top'] = False


def FHN_deriv(t, y):
    V,w = y 
    Vdot = (1/eps)*( V*(1-V)*(V-alpha) - w) + I
    wdot = V - w
    return [Vdot,wdot]


def get_V_W(t):   # For any t, find V and W on the stable periodic solution
    tt = np.remainder(t,T)
    dum = t_n-tt
    dum_arg = np.argmax(dum>0)
    tup = t_n[dum_arg]
    tdown = t_n[dum_arg-1]
    Vup = V_n[dum_arg]
    Vdown = V_n[dum_arg-1]
    V = Vdown + ((tt-tdown)/(tup-tdown))*(Vup-Vdown) # linear interpolation to find V(t)
    Wup = W_n[dum_arg]
    Wdown = W_n[dum_arg-1]
    W = Wdown + ((tt-tdown)/(tup-tdown))*(Wup-Wdown) # linear interpolation to find W(t)
    return V,W


def adj_deriv(t,y):
    V,W = get_V_W(t)        # For each t we need to find the value of V and W on the limit cycle
    y1dot = -(1/eps)*(-3*V*V + 2*V*(1+alpha) - alpha)*y[0] - y[1]
    y2dot = (1/eps)*y[0] + y[1]
    return [y1dot,y2dot]


alpha = 0.05
eps = 0.01  
I = 5

# First find the periodic solution
t0=0
tf=20
t = np.arange(t0,tf,0.001)
# Initial conditions
y0 = [0.1,0]
soln = solve_ivp(FHN_deriv, (t0, tf), y0, method='Radau',t_eval=t)
v,w = soln.y                                    # because I like giving things nicer names.
peaks,_ = find_peaks(v)
V_n = v[peaks[-2]:peaks[-1]]                    # Extract the solution over one period
W_n = w[peaks[-2]:peaks[-1]]
t_n = (t[peaks[-2]:peaks[-1]] - t[peaks[-2]])
T = t_n[-1]                                     # period of the solution 


# Now solve the adjoint equation backwards in time
t = np.linspace(0,-20,20000)
soln = solve_ivp(adj_deriv, (0, -20), y0, method='Radau',t_eval=t)
Va,Wa = soln.y                                  
peaks,_ = find_peaks(Wa)
Va = Va[peaks[-2]:peaks[-1]]                    # Extract the solution over one period
Wa = Wa[peaks[-2]:peaks[-1]]
ta = (t[peaks[-2]:peaks[-1]] - t[peaks[-2]])
Ta = -ta[-1]                                     # period of the solution

plt.figure(1)
plt.plot(ta/Ta,Va)
plt.show(1)

# Finally, normalise the adjoint solution
V_n = V_n[0:-1]
W_n = W_n[0:-1]
integrand = Va*((V_n*(1-V_n)*(V_n-alpha) - W_n)*(1/eps) + I) + Wa*(V_n-W_n)
norm = np.trapz(integrand,ta)/(T)
plt.figure(2)
plt.plot(ta,integrand)
plt.show(2)
print(norm)

np.savetxt('FHN_adjoint.dat',np.c_[ta/Ta,Va])
