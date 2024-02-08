
#   -------------------------------------------------------------------
# 
#    Use the method of lines to compute a traveling calcium wave. This
#    version uses the modal model of the IP3 receptor, but this can be
#    changed. Here, IP3 is also diffusing and has simple reaction terms.
# 
#    For Chapter 7, Section 7.6.1 of
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
    c = s[:N]
    ce = s[N:2*N]
    h = s[2*N:3*N]
    IP = s[3*N:]

    lam_c = Dc / (delx**2)
    lam_e = De / (delx**2)
    lam_p = DIP / (delx**2)

    out = coscrhs(t, s)
    
    Fc = np.zeros(N)
    Fce = np.zeros(N)
    Fh = np.zeros(N)
    FIP = np.zeros(N)
    
    Fc[0] =  lam_c * (-2*c[0] + 2*c[1]) + out[0]
    Fc[1:N-1] = lam_c * ( c[2:N] - 2*c[1:N-1] + c[0:N-2] ) + out[1:N-1]
    Fc[N-1] = lam_c*(-2*c[N-1] + 2*c[N-2]) + out[N-1]
    
    Fce[0] =  lam_e * (-2*ce[0] + 2*ce[1]) + out[N]
    Fce[1:N-1] = lam_e * ( ce[2:N] - 2*ce[1:N-1] + ce[0:N-2] ) + out[N+1:2*N-1]
    Fce[N-1] = lam_e*(-2*ce[N-1] + 2*ce[N-2]) + out[2*N-1]
    
    Fh = out[2*N:3*N]
    
    FIP[0] =  lam_p * (-2*IP[0] + 2*IP[1]) + out[3*N]
    FIP[1:N-1] = lam_p * ( IP[2:N] - 2*IP[1:N-1] + IP[0:N-2] ) + out[3*N+1:4*N-1]
    FIP[N-1] = lam_p*(-2*IP[N-1] + 2*IP[N-2]) + out[4*N-1]
    

    return np.concatenate((Fc, Fce, Fh, FIP))


# Define the right-hand sides of the ODEs
def coscrhs(t, s):
    c = s[:N]
    ce = s[N:2*N]
    h = s[2*N:3*N]
    IP = s[3*N:]

    phi_c = c**4 / (c**4 + Kc**4)
    phi_p = IP**2 / (Kp**2 + IP**2)
    phi_p_down = Kp**2 / (Kp**2 + IP**2)

    Jpm = Vpm * c**2 / (Kpm**2 + c**2)
    Jin = alpha0 + Vsocc * (Ksocc**4 / (Ksocc**4 + ce**4))
    h_inf = Kh**4 / (Kh**4 + c**4)
    tau = tau_max * Ktau**4 / (Ktau**4 + c**4)
    Jserca = Vserca * (c**2 - Kbar * ce**2) / (c**2 + kserca**2)

    beta = phi_p * phi_c * h
    alpha = phi_p_down * (1 - phi_c * h_inf)
    Po = beta / (beta + 0.4 * (beta + alpha))
    Jipr = kf * Po * (ce - c)

    dcdt = Jipr - Jserca + delta * (Jin - Jpm)
    dcedt = gamma * (Jserca - Jipr)
    dhdt = (h_inf - h) / tau
    dIPdt = Vplc - kp * IP

    return np.concatenate((dcdt, dcedt, dhdt, dIPdt))


tau_max = 100
delta = 1.5
Ktau = 0.1
tauP = 1
Kc = 0.2
Kh = 0.08
kf = 10
Vserca = 0.9
kserca = 0.2
Kbar = 0.00001957
Kp = 0.2
gamma = 5.5
Vpm = 0.11
Kpm = 0.3
alpha0 = 0.0027
Vsocc = 0.07
Ksocc = 8
Vplc = 0.02
kp = 1
Dc = 5  # Ca diffusion coefficient
De = 5  # ER Ca diffusion coefficient
DIP = 200  # IP3 diffusion coefficient

# Integration parameters
N = 300  # number of spatial grid points
L = 60  # length of domain
delx = L / N

# Initial condition
X = delx * np.arange(1, N + 1)
c0 = 0.08 * np.ones(N)
c0[:50] = 0.2  # initial calcium step on left side, to start a wave
ce0 = 15 * np.ones(N)
h0 = 0.45 * np.ones(N)
IP0 = (Vplc / kp) * np.ones(N)  # IP3: this is in the range where we expect a solitary pulse
init = np.concatenate((c0, ce0, h0, IP0))

# Time span
t_end = 10
n_out = 101
tspan = np.linspace(0, t_end, n_out)

# Integrate the system of ODEs
sol = odeint(pdeRHS, init, tspan)

# Plot two solutions at different times
plt.plot(X, sol[49,:N], label='t = 4.9', linewidth=2)
plt.plot(X, sol[69,:N], label='t = 6.9', linewidth=2)
plt.xlabel('x')
plt.ylabel(r'$\mathrm{Ca}^{++}(x,t)$')
plt.legend()
plt.show()

