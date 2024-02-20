
#   -------------------------------------------------------------------
# 
#    Calculate sensitivies in the five-compartment circulation model
# 
#    For Chapter 11, Section 11.5.4 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np

def getparams(V):
    # Define the nominal parameter values
    V0s = 0.94
    V0p = 0.26
    F = 80

    # Extract variables from the input vector V
    Psa, Ps, Psv, Ppa, Ppv, Vsa, Vsv, Vp, Q = V

    # Calculate parameters based on variables
    Rsa = (Psa - Ps) / Q
    Cld = Q / (F * Ppv)
    Csa = 2 * (Vsa - V0s) / (Psa + Ps)
    Rsv = (Ps - Psv) / Q
    Csv = 2 * Vsv / (Psv + Ps)
    Rp = (Ppa - Ppv) / Q
    Crd = Q / (Psv * F)
    Cp = 2 * (Vp - V0p) / (Ppa + Ppv)
    
    return [Csa, Csv, Cp, Rsa, Rsv, Rp, Cld, Crd]

def getpressures(P):
    # Extract parameters from the input vector P
    Csa, Csv, Cp, Rsa, Rsv, Rp, Cld, Crd, F = P

    # Define initial volumes
    V0s = 0.94
    V0p = 0.26
    Vt = 5
    Ve = Vt - V0s - V0p

    # Calculate alp
    alp = Rsv * (Csa + Csv / 2) + Rsa * Csa / 2 + Cp * Rp / 2

    # Calculate flow rate Q
    Q = Ve / (alp + Cp / (F * Cld) + (Csv + Csa) / (F * Crd))

    # Calculate pressures and volumes
    Psa = Q * (1 / (F * Crd) + Rsa + Rsv)
    Ps = Q * (1 / (F * Crd) + Rsv)
    Psv = Q / (F * Crd)
    Ppa = Q * (1 / (F * Cld) + Rp)
    Ppv = Q / (F * Cld)
    Vsa = V0s + Csa * (Psa + Ps) / 2
    Vsv = Csv * (Psv + Ps) / 2
    Vp = V0p + Cp * (Ppa + Ppv) / 2

    return np.array([Psa, Ps, Psv, Ppa, Ppv, Vsa, Vsv, Vp, Q])


# Specify the initial variables and set nominal values

V = np.array([100, 30, 2, 15, 5, 1, 3.5, 0.5, 5.6])

# Calculate nominal parameters
P = getparams(V)
P = np.append(P,80)  # Update the nominal rate

# Check the calculation
Vn = getpressures(P)

# Initialize sensitivity matrix
sens = np.zeros((9, 9))

# Now calculate sensitivities
for j in range(9):
    # Specify a perturbation for each parameter
    frac = 0.001
    del_p = P[j] * frac
    Pdel = P.copy()
    Pdel[j] = P[j] + del_p
    Vdel = getpressures(Pdel)
    sens[:, j] = (Vdel - V) / (frac * V)

# The sensitivity matrix
print(sens)


