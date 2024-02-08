
#   -------------------------------------------------------------------
# 
#    Use the method of lines to compute a traveling calcium wave. This
#    version uses a particularly simple wave model.
# 
#    For Chapter 7, Section 7.7.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Set default plotting parameters
plt.rcParams.update({
    'axes.labelsize': 20,
    'axes.linewidth': 2.0,
    'lines.linewidth': 2.0,
    'figure.figsize': (10, 8)
})

def pdeRHS(t, s):    
    # Extract variables
    c = s[:N]
    ce = s[N:2*N]

    # Calculate the reaction terms
    out = coscrhs(t, s)
    c_react = out[:N]
    ce_react = out[N:2*N]
    
    Fc = np.zeros(N)

    # Calculate the Ca diffusion terms.
    Fc[0] =  lam_c * (-2*c[0] + 2*c[1]) + c_react[0]
    Fc[1:N-1] = lam_c * ( c[2:N] - 2*c[1:N-1] + c[0:N-2] ) + c_react[1:N-1]
    Fc[N-1] = lam_c*(-2*c[N-1] + 2*c[N-2]) + c_react[N-1]
    
    Fce = ce_react   # no diffusion of Ca in the ER
    
    # Concatenate derivatives
    s_prime = np.concatenate((Fc, Fce))
    
    return s_prime

def coscrhs(t, s):
    # Extract variables
    c = s[:N]  # Calcium
    ce = s[N:2*N]  # ER calcium
    
    # Calculate open probability
    Po = p * (c**2 / (c**2 + ph1**2)) * (ph2 / (c + ph2))
    
    # Calculate fluxes
    Jipr = (kf * Po + a) * (ce - c)
    Jserca = ks * c
    Jpm = km * c
    Jin = a0 + a1 * p
    
    # Calculate derivatives
    Fc = Jipr - Jserca + delta * (Jin - Jpm)
    Fce = gam * (Jserca - Jipr)
    
    # Concatenate derivatives
    out = np.concatenate((Fc, Fce))
    
    return out


km = 1
ks = 20
a0 = 0.1
a1 = 0.1
ph1 = 2
ph2 = 1
a = 0.05
kf = 20
delta = 1
gam = 5
Dc = 1  # Ca diffusion coefficient
parlist = [2.5, 0.5]  # the bifurcation parameter
c0list = [0.5, 0.28]
ce0list = [5.2, 30]
axlist = [15, 35]

c = np.linspace(0.01, 1.5, 500)  # for the nullcline plots

########## Integrate the ODE for two different values of p ##########

for j in range(len(parlist)):
    p = parlist[j]  # the bifurcation parameter
    
    c0 = c0list[j]
    ce0 = ce0list[j]
    
    # make a phase portrait for each parameter on the list
    Po = p * (c**2 / (c**2 + ph1**2)) * (ph2 / (c + ph2))  # open probability
    
    # nullclines
    Jserca = ks * c
    Jpm = km * c
    Jin = a0 + a1 * p
    ce1 = (Jserca - delta * (Jin - Jpm)) / (kf * Po + a) + c
    ce2 = Jserca / (kf * Po + a) + c
    
    # integrate the ode
    init = [c0, ce0]  # initial data for the ode
    ct = c0 + ce0 / gam
    tstep = 0.05
    t_end = 100
    tspan = np.linspace(0, t_end, 2000)
    N = 1

    odesol = solve_ivp(coscrhs, [tspan[0], tspan[-1]], init, t_eval=tspan)
    csol = odesol.y[0,:]
    cesol = odesol.y[1,:]
    
    plt.figure(2*j - 1)  # a phase portrait
    plt.plot(c, ce1, '--', c, ce2, '--', csol, cesol, linewidth=2)
    plt.legend(['dc/dt=0', 'dc_e/dt=0'], bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel(r'$c$')
    plt.ylabel(r'$c_e$')
    plt.axis([0, 1.5, 0, axlist[j]])
    plt.title('p = {:.2f} µM'.format(p), fontsize=18)

    plt.show()
    
    

########## Integrate the PDE ##########

# Set parameters
N = 600            # Number of spatial grid points
L = 30             # Length of domain
h = L/N
lam_c = Dc / h ** 2
p = parlist[1]     # Bifurcation parameter

# Set initial data
X = h * np.arange(1, N+1)
V = np.exp(-3*X**2)  # Calcium

# Use the approximate steady state as initial data for ce
ce0 = 30.56 * np.ones(N)
init = np.concatenate((V, ce0))

# Define time points
tstep = 0.05  # Time between plots
t_end = 8     # Total time to run simulation
tspan = np.arange(0, t_end + tstep, tstep)

# Integrate the system of ODEs
S = solve_ivp(pdeRHS,  [tspan[0], tspan[-1]], init, t_eval=tspan)

# Plot a chosen solution
plt.figure(6)
plt.plot(X, S.y[:N,50], label='t=2.5', linewidth=2)
plt.plot(X, S.y[:N,100], label='t=5', linewidth=2)
plt.plot(X, S.y[:N,150],  label='t=7.5', linewidth=2)
plt.xlabel('x', fontsize=20)
plt.ylabel(r'$Ca^{2+}(x,t)$', fontsize=20)
plt.legend()

plt.figure(7)
plt.plot(X, S.y[N:2*N,50],X, S.y[N:2*N,100],X, S.y[N:2*N,150], linewidth=2)
plt.xlabel('x', fontsize=20)
plt.ylabel(r'$C_e(x,t)$', fontsize=20)

# Select an arbitrary point in the domain, and plot c against ce at that point
plt.figure(8)
plt.plot(c, ce1, '--', c, ce2, '--', S.y[10,:], S.y[N+10,:], linewidth=2)
plt.xlabel('Ca', fontsize=20)
plt.ylabel('C_e', fontsize=20)
plt.legend(['dc/dt=0', 'dc_e/dt=0'], loc='upper right')
plt.show()
