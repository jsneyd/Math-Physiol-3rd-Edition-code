#   -------------------------------------------------------------------
# 
#    Code for the Aplysia model.
# 
#    For Chapter 9, Exercise 9.10 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

from numpy import exp, arange
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ODE functions
def deRHS(sol, t):
    global rho, Kc, lam, cm, vna, vca, vk, vL, gna, gk, gkca, gL, C1, C2, taux, A, B, gca

    v, h, n, x, c = sol

    # Activation variables
    alpham = 0.1 * (50.0 - (C1 * v + C2)) / (exp((50.0 - (C1 * v + C2)) / 10.0) - 1.0)
    betam = 4.0 * exp((25.0 - (C1 * v + C2)) / 18.0)
    minf = alpham / (alpham + betam)

    alphah = 0.07 * exp((25.0 - (C1 * v + C2)) / 20.0)
    betah = 1.0 / (exp((55.0 - (C1 * v + C2)) / 10.0) + 1.0)
    hinf = alphah / (alphah + betah)
    tauh = 1.0 / (alphah + betah)

    alphan = 0.01 * (55.0 - (C1 * v + C2)) / (exp((55.0 - (C1 * v + C2)) / 10.0) - 1.0)
    betan = 0.125 * exp((45.0 - (C1 * v + C2)) / 80.0)
    ninf = alphan / (alphan + betan)
    taun = 1.0 / (alphan + betan)

    xinf = 1.0 / (exp(A * (B - v)) + 1.0)

    # Ionic currents
    ina = gna * (minf**3) * h * (v - vna)
    ica = gca * x * (v - vca)
    iL = gL * (v - vL)
    ik = (gk * (n**4) + gkca * c / (0.5 + c)) * (v - vk)

    # Differential Equations
    vp = (-ina - ica - ik - iL) / cm
    hp = lam * (hinf - h) / tauh
    np = lam * (ninf - n) / taun
    xp = (xinf - x) / taux
    cp = rho * (Kc * x * (vca - v) - c)

    return [vp, hp, np, xp, cp]

# Global parameters
rho = 0.0003
Kc = 0.0085
lam = 0.08
cm = 1.0
vna = 30.0
vca = 140.0
vk = -75.0
vL = -40.0
gna = 4.0
gca = 0.004
gk = 0.3
gkca = 0.03
gL = 0.003
C1 = 1.209
C2 = 78.714
taux = 235.0
A = 0.15
B = -50.0

# Initial conditions
v0 = -70
h0 = 0.012925
n0 = 0.61732
x0 = 0.5
c0 = 0.6
u0 = [v0, h0, n0, x0, c0]

total = 40000

tstep = 1
# Specify the output points
tspan = arange(0, total + tstep, tstep)

# Solve the ODE
sol = odeint(deRHS, u0, tspan)

plt.figure(1)
plt.plot(tspan / 1000, sol[:, 0])
plt.ylabel('V (mV)')
plt.xlabel('t (s)')

plt.show()



