# Use shooting to find the periodic solution of the adjoint problem for a FHN oscillator

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
plt.rcParams['axes.spines.right'] = False  # Me being pedantic and removing the box around the graph
plt.rcParams['axes.spines.top'] = False

def FHN_deriv(t, y):
    v,w = y 
    vdot = (1/eps)*( v*(1-v)*(v-alpha) - w) + I
    wdot = v - w
    return [vdot,wdot]

def get_v_w(t):   # For any t, find V and W on the stable periodic solution
    tt = np.remainder(t,T)
    dum = t_n-tt
    dum_arg = np.argmax(dum>0)
    tup = t_n[dum_arg]
    tdown = t_n[dum_arg-1]
    vup = v_n[dum_arg]
    vdown = v_n[dum_arg-1]
    v = vdown + ((tt-tdown)/(tup-tdown))*(vup-vdown) # linear interpolation to find V(t)
    wup = w_n[dum_arg]
    wdown = w_n[dum_arg-1]
    w= wdown + ((tt-tdown)/(tup-tdown))*(wup-wdown) # linear interpolation to find W(t)
    return v,w

def adj_deriv(t,y):
    v,w = get_v_w(t)        # For each t we need to find the value of v and w on the limit cycle
    y1dot = -(1/eps)*(-3*v*v + 2*v*(1+alpha) - alpha)*y[0] - y[1]
    y2dot = (1/eps)*y[0] + y[1]
    return [y1dot,y2dot]

def FHN_objective(y):
    vinit,winit,T = y
    t = np.linspace(0,T,2000)
    y0 = [vinit,winit]
    soln = solve_ivp(FHN_deriv, (0, T), y0, method='Radau',t_eval=t,rtol=1e-8,atol=1e-8)
    v,w = soln.y
    return (v[1]-v[-1])**2 + (w[1]-w[-1])**2

def adjoint_objective(Y):
    Vinit,Winit = Y
    t = np.linspace(T,0,2000)
    Y0 = [Vinit,Winit]
    soln = solve_ivp(adj_deriv, (T, 0), Y0, method='Radau',t_eval=t,rtol=1e-8,atol=1e-8)
    V,W = soln.y
    return (V[1]-V[-1])**2 + (W[1]-W[-1])**2


alpha = 0.05
eps = 0.01  
I = 5

###############################################

# Find the periodic solution of FHN by shooting.

y0 = [0.1,0,1.2]  # in the order v0, w0, T. Arbitrary choices, but they seem to work
res = minimize(FHN_objective, y0)
v0,w0,T = res.x

# Now construct the FHN periodic solution, to be used later
t = np.linspace(0,T,2000)
y0 = [v0,w0]
soln = solve_ivp(FHN_deriv, (0, T), y0, method='Radau',t_eval=t,rtol=1e-8,atol=1e-8)
v_n,w_n = soln.y   # Here's our periodic solution
t_n = soln.t
#plt.plot(soln.t,v_n)

###############################################

# get the initial condition for the shooting method

t = np.linspace(3*T,0,2000)
y0 = [v0,w0]
soln = solve_ivp(adj_deriv, (3*T, 0), y0, method='Radau',t_eval=t,rtol=1e-8,atol=1e-8)
#plt.plot(soln.t,soln.y[0])
Va,Wa = soln.y

# Find the periodic solution of the adjoint by shooting.
Y0 = [Va[-1],Wa[-1]]
res = minimize(adjoint_objective, Y0)
V0,W0 = res.x

# construct and plot the adjoint solution
t = np.linspace(T,0,2000)
Y0 = [V0,W0]
soln = solve_ivp(adj_deriv, (T, 0), Y0, method='Radau',t_eval=t,rtol=1e-8,atol=1e-8)
Va,Wa = soln.y
#plt.plot(t_n,v_n,soln.t,Va)
Va = np.flip(Va)            # so it is the proper time order
Wa = np.flip(Wa)

###############################################

# Finally, normalise the adjoint solution

integrand = Va*((v_n*(1-v_n)*(v_n-alpha) - w_n)*(1/eps) + I) + Wa*(v_n-w_n)
norm = np.trapz(integrand,t_n)/T
plt.figure(2)
plt.plot(t_n,integrand)
plt.xlabel('time')
plt.ylabel('normalisation integrand')
plt.show(2)
print(norm)


# Plot the normalised adjoint solution
plt.figure(3)
plt.plot(t_n/T,Va/norm)

np.savetxt('FHN_adjoint.dat',np.c_[t_n/T,Va/norm])
