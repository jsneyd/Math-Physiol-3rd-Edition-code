#  -------------------------------------------------------------------
#
#   Solve the model of bursting in GnRH neurons.
#   Simplified version of the model of Duan et al, 
#           Journal of Theoretical Biology 276 (2011) 22–34
#
#   For Chapter 9, Section 9.2 of
#   Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#   Written by James Keener and James Sneyd.
# 
#  ------------------------------------------------------------------- 


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Function defining the model
def model(t,S):
    V, Hnaf, Nkm, c, ce, h, Ow, Ow_star = S
    
    Mnaf_inf = 1/(1+np.exp(-(V+40.0)/4.3))
    Hnaf_inf = 1/(1+np.exp((V+66.1)/10.8))
    Nkm_inf = 1/(1+np.exp(-(V+37)/4))
    Mcal_inf = 1/(1+np.exp(-(V+30)/2))
    
    Thnaf = 75/(np.exp((V+80)/19)+2*np.exp(-2*(V+80)/19))
    Tnkm = 11.5/(np.exp((V+30)/15)+np.exp(-(V+30)/15))
    
    Inaf = gnaf*Mnaf_inf**3*Hnaf*(V-Vna)
    Ical = gcal*Mcal_inf**2*(V-Vca)
    Ileak = gleak*(V-Vleak)
    
    Ikm = gkm*Nkm*(V-Vk)
    Isk = gsk*(c**(n_sk)/(c**(n_sk)+k_sk**(n_sk)))*(V-Vk)
    Iw = gw*(Ow+Ow_star)*(V-Vk)
    
    Iionic = Inaf+Ikm+Ical+Ileak
    
    # calcium submodel
    phi_c=c**4/(c**4+Kc**4)
    phi_p=IP3**2/(Kp**2+IP3**2)
    phi_p_down=Kp**2/(Kp**2+IP3**2)
    
    h_inf=Kh**4/(Kh**4+c**4)
    tau = tau_max*Ktau**4/(Ktau**4+c**4)
    beta = phi_p*phi_c*h
    alpha = phi_p_down*(1-phi_c*h_inf)
    Po=beta/(beta + 0.4*(beta+alpha))
    
    JIPR= Kf*Po*(ce-c);
    Jserca = Vserca*(c*c-Kbar*ce*ce)/(c*c+kserca*kserca)
    Jin = -alpha0*(Ical)+alpha1*IP3
    Jpm = Vp*c**2/(Kp**2+c**2)
    
    dVdt = (-1/Cm)*(Iionic+Isk+Iw)
    dHnafdt = (Hnaf_inf-Hnaf)/Thnaf
    dNkmdt = (Nkm_inf-Nkm)/Tnkm
    dcdt = delta*(Jin - Jpm) + (JIPR - Jserca);
    dhdt  = 0.001*(h_inf-h)/tau   # Note the factor of 1/1000, since time is in ms, not s.
    dcedt = gamma*(Jserca - JIPR)
    dOwdt = -Ow*k_11-Ow*k22+k11*c*(1-Ow-Ow_star)
    dOw_stardt = k22*Ow-k33*Ow_star
        
    return [
        dVdt, dHnafdt, dNkmdt,
        dcdt, dcedt, dhdt,
        dOwdt, dOw_stardt
    ]



# Parameters
IP3 = 0.15
Kf = 1
tau_max = 10
Ktau = 0.1
tauP = 0.5
Kc = 0.2
Kh = 0.08
Kp = 0.425
delta = 0.5
gamma = 27
alpha0 = 4.8e-3
alpha1 = 2e-5
Vp = 0.0042
Vserca = 0.1
kserca = 0.2
Kbar = 0.00001957
gnaf = 280
gsk = 0.3
gw = 950
gcal = 0.05
gkm = 8
gleak = 0.04
Cm = 16
Vna = 60
Vca = 100
Vk = -80
Vleak = 100
k_sk = 1
n_sk = 2
k11 = 1e-7
k_11 = 1.2
k22 = 0.5
k33 = 3e-5

# Initial conditions
init = [
    -47.439744033440896,
    0.150875059704950,
    0.068499547036174,
    0.042896908631945,
    10,
    0.061113223260552,
    2.523211e-9,
    0.000049694789395
]

# Time points
t = np.linspace(0, 60000, 100000)


# Solve the ODEs
sol = solve_ivp(model, [0,60000], init, method='Radau',t_eval=t,max_step=1,rtol=1e-10,atol=1e-10)

# Extract variables
V = sol.y[0, :]
c = sol.y[3,:]
ce = sol.y[4,:]

# Plotting
plt.figure()
plt.plot(t/1000, V, label='V')
plt.xlabel('Time')
plt.ylabel('V (mV)')
plt.axis([52, 58, -65, 30])  # Set the x-axis and y-axis limits for the V plot


plt.twinx()
# Plot c on the secondary y-axis
plt.plot(t / 1000, c, label=r'$[{\rm Ca}^{2+}]$', color='red')
plt.ylabel(r'$[{\rm Ca}^{2+}]$ ($\mu$M)')
plt.legend()

