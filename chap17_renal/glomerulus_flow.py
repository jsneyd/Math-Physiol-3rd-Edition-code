import numpy as np
import matplotlib.pyplot as plt

# Global variables (using a dictionary to simulate MATLAB's global variables)
params = {}

# Define constants
P1 = 60
P2 = 18
Pa = 100
Pe = 18
Pd = 18  # Same as P2, so that Rd = 0
pii = 25  # mm Hg
Qi = 650
Qd = 125
Qe0 = Qi - Qd
alpha = pii / (P1 - P2)

# Calculate Ra, Re, Rd, and KfL
Ra = (Pa - P1) / Qi
Re = (P1 - Pe) / Qe0
Rd = (P2 - Pd) / Qd
Qe = Qe0
KfL = -(Qe / Qi + alpha * np.log((Qe / Qi - alpha) / (1 - alpha)) - 1) * (alpha * Qi) / pii

# Store parameters in the dictionary
params['KfL'] = KfL
params['pii'] = pii
params['Ra'] = Ra
params['Re'] = Re
params['Qi'] = Qi
params['Pe'] = Pe
params['P2'] = P2
params['Pa'] = Pa

# List of Pa values to iterate over
Palist = np.arange(20, 161)

# Store results for plotting
Qe_list = []

# Define the function for F(Qe)
def F(Qe, params):
    Qi = (params['Pa'] - params['Pe'] - params['Re'] * Qe) / params['Ra']
    P1 = -params['Ra'] * Qi + params['Pa']
    alpha = params['pii'] / (P1 - params['P2'])
    out = Qe + alpha * Qi * np.log((Qe / Qi - alpha) / (1 - alpha)) - Qi + params['KfL'] * params['pii'] / alpha
    return out

# Bisection method for finding the root
def bisect(f, a, b, params, tol=1e-6, max_iter=25):
    ul, uu = a, b
    fl, fu = f(ul, params), f(uu, params)

    for _ in range(max_iter):
        u = (ul + uu) / 2
        fc = f(u, params)

        if np.abs(fc) < tol:
            return u
        
        if fc * fl >= 0:
            ul, fl = u, fc
        else:
            uu, fu = u, fc
    
    return u

# Calculate Qe for each Pa in Palist
for Pa in Palist:
    params['Pa'] = Pa
    qemx = (Pa - Pe) / (Ra + Re)
    qemin = (P2 * Ra - Pe * Ra - Re * pii + np.sqrt((P2**2 - 2 * P2 * Pe + Pe**2) * Ra**2
               - 2 * (P2 - 2 * Pa + Pe) * Re * Ra * pii + Re**2 * pii**2)) / (2 * Re * Ra) + 0.1

    Qe_value = bisect(F, qemin, qemx, params)
    Qe_list.append(Qe_value)

# Convert lists to numpy arrays for plotting
Qe_array = np.array(Qe_list)
Qi_array = (Palist - Pe - Re * Qe_array) / Ra
Qd_array = Qi_array - Qe_array

# Filter Qd > 0
valid_indices = Qd_array > 0

# Plotting
plt.figure(1)
plt.plot(Palist[valid_indices], Qd_array[valid_indices], label='Q_d (ml/min)')
plt.xlabel('Arterial pressure, P_a (mm Hg)')
plt.ylabel('Glomerular filtration rate, Q_d (ml/min)')
plt.legend(loc="upper left")

plt.twinx()
plt.plot(Palist[valid_indices], Qe_array[valid_indices], color='orange', label='Q_e (ml/min)')
plt.ylabel('Efferent arteriole flow rate, Q_e (ml/min)')
plt.legend(loc="upper right")

plt.show()
