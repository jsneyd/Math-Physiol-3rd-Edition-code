
#   -------------------------------------------------------------------
# 
#    Use the method of lines to compute an intercellular calcium wave. This
#    version uses the modal model of the IP3 receptor, but this can be
#    changed. Here, IP3 is also diffusing and has simple reaction terms.
# 
#    This version does not have jump conditions at the intercellular boundaries, 
#    as Python doesn't seem to like DAEs. So the jump conditions have been 
#    converted to additional ODEs instead.
# 
#    For Chapter 7, Section 7.10.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# The right-hand side for PDE (MoL) simulation:
def pde_rhs(t, s):
    # There are four variables:  c, ce, h, IP
    # only c, ce, and IP are diffusing

    c = s[: N]  # calcium
    ce = s[N : 2 * N]  # ER calcium
    h = s[2 * N : 3 * N]  # h
    IP = s[3 * N :]  # IP3

    lam_c = Dc / delx ** 2
    lam_e = De / delx ** 2
    lam_p = DIP / delx ** 2

    out = coscrhs(t, s)
    sc = np.concatenate([np.array([1]), 2 * np.ones(N - 2), np.array([1])])

    Fc = (
        lam_c * (-sc * c + np.concatenate([np.array([0]), c[:-1]]) + np.concatenate([c[1:], np.array([0])]))
        + out[: N]
    )
    Fce = (
        lam_e
        * (-sc * ce + np.concatenate([np.array([0]), ce[:-1]]) + np.concatenate([ce[1:], np.array([0])]))
        + out[N : 2 * N]
    )
    Fh = out[2 * N : 3 * N]  # reaction only, no diffusion
    FIP = (
        lam_p
        * (-sc * IP + np.concatenate([np.array([0]), IP[:-1]]) + np.concatenate([IP[1:], np.array([0])]))
        + out[3 * N :]
    )

    # now put in the intercellular boundaries. Differentiate the jump conditions
    # to get these ODEs for the boundary points.
    for j in range(Ncell - 1):
        d = int((j + 1) * nc)

        Fc[d] = (Fc[d-1]*(1+muc) + muc*Fc[d+2])/(2*muc+1)
        Fc[d + 1] = (Fc[d+2]*(1+muc) + muc*Fc[d-1])/(2*muc+1)
        
        Fce[d] = (Fce[d-1]*(1+mue) + mue*Fce[d+2])/(2*mue+1)
        Fce[d + 1] = (Fce[d+2]*(1+mue) + mue*Fce[d-1])/(2*mue+1)
        
        FIP[d] = (FIP[d-1]*(1+muIP) + muIP*FIP[d+2])/(2*muIP+1)
        FIP[d + 1] = (FIP[d+2]*(1+muIP) + muIP*FIP[d-1])/(2*muIP+1)

    return np.concatenate([Fc, Fce, Fh, FIP])


# The right-hand sides of the ODEs
def coscrhs(t, s):
    c = s[: N]  # calcium
    ce = s[N : 2 * N]  # ER calcium
    h = s[2 * N : 3 * N]  # h
    IP = s[3 * N :]  # IP3

    # the equations for the modal model
    phi_c = c ** 4 / (c ** 4 + Kc ** 4)
    phi_p = IP ** 2 / (Kp ** 2 + IP ** 2)
    phi_p_down = Kp ** 2 / (Kp ** 2 + IP ** 2)

    Jpm = Vpm * c ** 2 / (Kpm ** 2 + c ** 2)
    Jin = alpha0 + Vsocc * (Ksocc ** 4 / (Ksocc ** 4 + ce ** 4))
    h_inf = Kh ** 4 / (Kh ** 4 + c ** 4)
    tau = tau_max * Ktau ** 4 / (Ktau ** 4 + c ** 4)
    Jserca = Vserca * (c * c - Kbar * ce * ce) / (c * c + kserca * kserca)

    beta = phi_p * phi_c * h
    alpha = phi_p_down * (1 - phi_c * h_inf)
    Po = beta / (beta + 0.4 * (beta + alpha))
    Jipr = kf * Po * (ce - c)

    dcdt = Jipr - Jserca + delta * (Jin - Jpm)
    dcedt = gamma * (Jserca - Jipr)
    dhdt = (h_inf - h) / tau
    dIPdt = Vplc * c ** 2 / (Kplc ** 2 + c ** 2) - kp * IP

    return np.concatenate([dcdt, dcedt, dhdt, dIPdt])



# Set parameters
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

Vplc = 0.02  # for a regenerating wave
#Vplc = 0    # for no regeneration
Kplc = 0.1
kp = 0.1

Dc = 5  # Ca diffusion coefficient
De = 5  # ER Ca diffusion coefficient
DIP = 200  # IP3 diffusion coefficient
FIP = 10  # intercellular IP3 permeability
Fc = 0  # intercellular Ca permeability
Fe = 0  # intercellular ER Ca permeability

# Geometry parameters
Ncell = 5  # number of cells
N = 500  # number of spatial grid points. Must be divisible by Ncell.
nc = N/Ncell  # number of spatial grid points in each cell.
L = 150  # length of domain
delx = L / N

muc = Fc * delx / Dc  # used in the cell boundary terms
mue = Fe * delx / De
muIP = FIP * delx / DIP

# Set initial data
# all variables start flat except for the IP3 pulse at the left
X = delx * np.arange(1, N + 1)
c0 = 0.08 * np.ones(N)
ce0 = 15 * np.ones(N)
h0 = 0.45 * np.ones(N)
IP0 = np.zeros(N)
IP0[:25] = 5  # IP3 stimulus at left of domain
init = np.concatenate([c0, ce0, h0, IP0])

# Specify the output points
t_end = 20  # total time to run simulation
n_out = 501  # number of output times
tspan = np.linspace(0, t_end, n_out)

# Integrate the system of ODEs
sol = solve_ivp(
    pde_rhs,
    t_span=(tspan[0], tspan[-1]),
    y0=init,
    method="BDF",
    t_eval=tspan,
    jac=None,
    rtol=1e-6,
    atol=1e-8
)

c = sol.y[:N]
ce = sol.y[N : 2 * N]
IP = sol.y[3 * N :]

# Plot a movie
plt.figure(1)
for i in range(n_out):
    plt.plot(X, c[:, i], linewidth=2)
    plt.ylim([0, 1])
    plt.xlim([0, L])
    plt.xlabel("x", fontsize=20)
    plt.ylabel(r"$[Ca^{2+}] (\mu {\rm M})$", fontsize=20)
    plt.pause(0.1)


plt.figure(2)
plt.imshow(np.transpose(c), aspect='auto', extent=[0, L, t_end, 0], cmap='viridis')
plt.colorbar(label=r'$Ca^{2+}$ concentration')
plt.xlabel("x ($\mu m$)")
plt.ylabel("time (sec)")
plt.show()





