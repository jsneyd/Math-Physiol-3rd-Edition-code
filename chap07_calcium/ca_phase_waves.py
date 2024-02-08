
#   -------------------------------------------------------------------
# 
#    Use the method of lines to compute a traveling calcium wave. 
#    Here, IP3 is also diffusing and has simple reaction terms,
#    so this version is used to demonstrate phase waves.
# 
#    For Chapter 7, Section 7.9 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Set default plotting parameters
plt.rcParams.update({
    'axes.labelsize': 20,
    'axes.linewidth': 2.0,
    'lines.linewidth': 2.0,
    'figure.figsize': (10, 8)
})

# Define the right-hand side of the PDE
def pdeRHS(s,t):
    IP = s[:N]
    c = s[N:2*N]
    y = s[2*N:3*N]

    lam_c = Dc / (delx**2)
    lam_p = DIP / (delx**2)

    out = coscrhs(s,t)
    
    FIP = np.zeros(N)
    Fc = np.zeros(N)
    Fy = np.zeros(N)
    
    FIP[0] =  lam_p * (-2*IP[0] + 2*IP[1]) + out[0]
    FIP[1:N-1] = lam_p * ( IP[2:N] - 2*IP[1:N-1] + IP[0:N-2] ) + out[1:N-1]
    FIP[N-1] = lam_p*(-2*IP[N-1] + 2*IP[N-2]) + out[N-1]
    
    Fc[0] =  lam_c * (-2*c[0] + 2*c[1]) + out[N]
    Fc[1:N-1] = lam_c * ( c[2:N] - 2*c[1:N-1] + c[0:N-2] ) + out[N+1:2*N-1]
    Fc[N-1] = lam_c*(-2*c[N-1] + 2*c[N-2]) + out[2*N-1]
    
    Fy = out[2*N:3*N]
    
    return np.concatenate((FIP, Fc, Fy))


# Define the right-hand sides of the ODEs
def coscrhs(s,t):
    # Evaluate the ODE dynamics
    IP = s[:N]        # IP3
    c = s[N:2*N]  # Calcium
    y = s[2*N:3*N]  # Inactivation
    ce = gam*(ct-c) # ER calcium

    Po = (IP * c * (1 - y) / ((IP + K1) * (c + K5)))**3  # Open probability
    ph1 = (km4 * K1 * K2 + km2 * K4 * IP) * c / (K4 * K2 * (IP + K1))
    ph2 = (km2 * IP + km4 * K3) / (K3 + IP)
    Jipr = kf * Po * (ce - c)
    Jserca = Vserca * (c**2 - Kbar * ce**2) / (Kserca**2 + c**2)

    Fp = np.zeros(N)  # IP3 is not reacting
    Fc = Jipr - Jserca
    Fy = ph1 * (1 - y) - ph2 * y

    return np.concatenate((Fp, Fc, Fy))



# parameters
Vserca = 0.9
Kbar = 2.e-7
kf=1.11
Kserca = 0.1
 
gam = 5.5
p = 0      # background ip3  
ct = 1.85  # total calcium
k1 = 400
k2 = 0.2
k3 = 400
k4 = 0.2
k5 = 20
km1 = 52
km2 = 0.21
km3 = 377.2
km4 = 0.0289
km5 = 1.64
K1=km1/k1
K5=km5/k5
K2 = km2/k2
K3=km3/k3
K4=km4/k4

Dc = 0.001  # Ca diffusion coefficient
DIP = 0.05  # IP3 diffusion coefficient

####  first integrate the odes ####

# variables are P, C, and y, in that order. Ce is found by conservation.
# Set initial conditions
N = 1   # just one space point to solve the odes
init = [0.35, 0.2, 1.0]  # Initial values for P, C, and y

# Set time span
tstep = 1
t_end = 200
tspan = np.arange(0, t_end + tstep, tstep)

# Integrate the system of ODEs
S = odeint(coscrhs, init, tspan)

# Extract individual variables
P = S[:, 0]
C = S[:, 1]
y = S[:, 2]
Ce = gam * (ct - C)

# # Plot the results
# plt.plot(tspan, P, label='IP3')
# plt.plot(tspan, C, label='Ca++')
# plt.plot(tspan, Ce, label='Ce')
# plt.plot(tspan, y, label='y')
# plt.xlabel('Time')
# plt.ylabel('Concentration')
# plt.legend()
# plt.show()



####  next integrate the pdes ####

# Integration parameters
N = 100  # number of spatial grid points
L = 5  # length of domain
delx = L / N

# Initial condition
X = delx * np.arange(1, N + 1)
c0 = C[-1] * np.ones(N)
y0 = y[-1] * np.ones(N)
IP0 = P[-1]*np.ones(N) + 3*np.exp(-3*X*X);  # initial ip3 concentration has a blob at one end
init = np.concatenate((IP0, c0, y0))

# Time span
tstep = 0.5  # time between plots
t_end = 150  # total time to run simulation
tspan = np.arange(0, t_end + tstep, tstep)

# Integrate the system of ODEs
sol = odeint(pdeRHS, init, tspan)


# If you want to make a simple movie of the wave, uncomment the next block

# # Loop through each time step and plot to make a simple movie
# for j in range(len(tspan)):
#     # Create a new figure
#     plt.figure(3)
    
#     # Plot Ca++(x, t)
#   #  plt.subplot(2, 1, 1)
#     plt.plot(X, sol[j, N:2*N], linewidth=2)
#     plt.xlabel('x', fontsize=20)
#     plt.ylabel('Ca++(x,t)', fontsize=20)
#     plt.axis([0, max(X), 0, 0.5])  # Set the axis limits
#     plt.show()



# Create a contour plot
plt.figure(7)
contour = plt.contour(X, tspan, sol[:, N:2*N], linewidths=2)
plt.xlabel('x')
plt.ylabel('t')
plt.title('Contour Plot of Ca++')
plt.clabel(contour, inline=True, fontsize=8)
plt.show()

    

